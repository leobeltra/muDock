#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

#include <oneapi/tbb/parallel_pipeline.h>
#include <oneapi/tbb/global_control.h>
#include <memory>
#include <vector>
#include <iostream>

#include <mpi.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>

#include <cstdint>

namespace mudock {

    static inline void print_ligand(const static_molecule& ligand) {
        std::cout << ligand.properties.get(property_type::NAME) << " "
                  << ligand.properties.get(property_type::SCORE) << "\n";
    }

    void run_tbb_pipeline(std::istream& in,
         const std::vector<std::string>& configurations,
         const knobs& knobs,
         genetic_adt_pipeline& pipeline, 
         std::size_t end,
         std::size_t max_tokens)   
    {
        using mol_vec = parser_filter::mol_vec;

        auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
        auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
        input_queue->initialize(1000);
        output_queue->initialize(1000);

        std::thread writer([&]{
            bool ok=false;
            std::string buf;
            buf.reserve(1<<20);

            std::size_t lines = 0;
            for (auto x = output_queue->dequeue(ok); ok; x = output_queue->dequeue(ok)) {
                buf += x->properties.get(property_type::NAME);
                buf += ' ';
                buf += x->properties.get(property_type::SCORE);
                buf += '\n';

                if (++lines % 4096 == 0) {
                    std::cout << buf;
                    buf.clear();
                }
            }
            if (!buf.empty()) std::cout << buf;
        });

        // auto drain_ready_nb = [&]() {
        //     mudock::safe_queue<mudock::static_molecule>::value_ptr_type x;
        //     while (output_queue->try_dequeue(x)) {     
        //         print_ligand(*x);
        //     }
        // };

        // auto drain_ready = [&]() -> std::size_t {
        //     std::size_t counter = 0;
        //     bool is_retrieved = false;
        //     while (auto x = output_queue->dequeue(is_retrieved)) {
        //         print_ligand(*x);
        //         ++counter;
        //     }
            
        //     return counter;
        // };

        threadpool pool;
        manager(configurations, pool, knobs, input_queue, output_queue, pipeline);
        info("Manager done: workers created");

        const std::size_t effective_max_tokens = max_tokens == 0 ? 1 : max_tokens;

        oneapi::tbb::parallel_pipeline(
            effective_max_tokens,
            oneapi::tbb::make_filter<void, std::string>(
            oneapi::tbb::filter_mode::serial_in_order,
            stream_filter(in, end)
            )
            &
            oneapi::tbb::make_filter<std::string, mol_vec>(
            oneapi::tbb::filter_mode::parallel,
            parser_filter()
            )
            &
            oneapi::tbb::make_filter<mol_vec, void>(
                oneapi::tbb::filter_mode::serial_out_of_order,
            [&](mol_vec molecules) {
                bool is_produced = false;
                for (auto& p : molecules) {
                    if (p) {
                        input_queue->enqueue(p, is_produced);
                    }
                    if (!is_produced) break;
                }
            }
            )
            // &
            // oneapi::tbb::make_filter<std::size_t, void>(
            //         oneapi::tbb::filter_mode::parallel,
            //     [&] (std::size_t) {
            //         drain_ready_nb();
            //     }
            // )
        );

        info("Pipeline done: closing input queue");
        input_queue->send_terminate_signal();
        pool.wait();
        output_queue->send_terminate_signal();
        writer.join();

        info("Output drained to file");
    }

} // namespace mudock

// version with tbb writer stage
        
// #include <mudock/tbb_implementation/stream_filter.hpp>
// #include <mudock/tbb_implementation/parser_filter.hpp>
// #include <mudock/tbb_implementation/tbb_pipeline.hpp>
// #include <mudock/compute/manager.hpp>
// #include <mudock/compute/safe_queue.hpp>
// #include <mudock/compute/safe_stack.hpp>
// #include <mudock/molecule.hpp>
// #include <mudock/mudock.hpp>

// #include <oneapi/tbb/parallel_pipeline.h>
// #include <memory>
// #include <vector>
// #include <iostream>

// namespace mudock {

//     static inline void print_ligand(const static_molecule& ligand) {
//         std::cout << ligand.properties.get(property_type::NAME) << " "
//                   << ligand.properties.get(property_type::SCORE) << "\n";
//     }

//     void run_tbb_pipeline(std::istream& in,
//                           const std::vector<std::string>& configurations,
//                           const knobs& knobs,
//                           genetic_adt_pipeline& pipeline,
//                           std::size_t end,
//                           std::size_t max_tokens)
//     {
//         using mol_vec = parser_filter::mol_vec;

//         auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
//         auto output_queue = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();

//         input_queue->initialize(1000);

//         auto drain_ready_nb = [&]() -> std::size_t {
//             std::size_t drained = 0;
//             std::string buf;
//             buf.reserve(1 << 20);

//             while (auto x = output_queue->dequeue()) {
//                 buf += x->properties.get(property_type::NAME);
//                 buf += ' ';
//                 buf += x->properties.get(property_type::SCORE);
//                 buf += '\n';
//                 ++drained;

//                 if ((drained % 4096) == 0) {
//                     std::cout << buf;
//                     buf.clear();
//                 }
//             }

//             if (!buf.empty()) {
//                 std::cout << buf;
//             }

//             return drained;
//         };

//         auto drain_all_blocking = [&]() {
//             std::string buf;
//             buf.reserve(1 << 20);

//             std::size_t lines = 0;
//             for (auto x = output_queue->dequeue_wait(); x; x = output_queue->dequeue_wait()) {
//                 buf += x->properties.get(property_type::NAME);
//                 buf += ' ';
//                 buf += x->properties.get(property_type::SCORE);
//                 buf += '\n';

//                 if (++lines % 4096 == 0) {
//                     std::cout << buf;
//                     buf.clear();
//                 }
//             }

//             if (!buf.empty()) {
//                 std::cout << buf;
//             }
//         };

//         threadpool pool;
//         manager(configurations, pool, knobs, input_queue, output_queue, pipeline);
//         info("Manager done: workers created");

//         const std::size_t effective_max_tokens = max_tokens == 0 ? 1 : max_tokens;

//         oneapi::tbb::parallel_pipeline(
//             effective_max_tokens,
//             oneapi::tbb::make_filter<void, std::string>(
//                 oneapi::tbb::filter_mode::serial_in_order,
//                 stream_filter(in, end)
//             )
//             &
//             oneapi::tbb::make_filter<std::string, mol_vec>(
//                 oneapi::tbb::filter_mode::parallel,
//                 parser_filter()
//             )
//             &
//             oneapi::tbb::make_filter<mol_vec, std::size_t>(
//                 oneapi::tbb::filter_mode::serial_out_of_order,
//                 [&](mol_vec molecules) -> std::size_t {
//                     bool is_produced = false;
//                     std::size_t enq = 0;

//                     for (auto& p : molecules) {
//                         if (!p) continue;

//                         input_queue->enqueue(p, is_produced);
//                         if (!is_produced) break;
//                         ++enq;
//                     }

//                     return enq;
//                 }
//             )
//             &
//             oneapi::tbb::make_filter<std::size_t, void>(
//                 oneapi::tbb::filter_mode::serial_out_of_order,
//                 [&](std::size_t) {
//                     // drain data
//                     drain_ready_nb();
//                 }
//             )
//         );

//         info("Pipeline done: closing input queue");
//         input_queue->send_terminate_signal();

//         pool.wait();

//         output_queue->close();

//         drain_all_blocking();

//         info("Output drained to file");
//     }

// } // namespace mudock