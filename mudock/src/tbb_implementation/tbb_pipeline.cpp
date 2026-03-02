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
    // =========================
    // consumer thread (NEW)
    // // =========================
    // std::atomic<bool> done{false};

    // std::thread consumer([&]{
    //     std::ofstream out("results/rank_" + std::to_string(rank) + ".csv");
    //     while (!done.load()) {
    //         auto x = output_queue->dequeue();
    //         if (x) {
    //             out << x->properties.get(property_type::NAME) << ","
    //                 << x->properties.get(property_type::SCORE) << "\n";
    //         } else {
    //             std::this_thread::sleep_for(std::chrono::microseconds(50));
    //         }
    //     }
    // });

        auto drain_ready = [&]() -> std::size_t {
        std::size_t counter = 0;
        while (auto x = output_queue->dequeue()) {
            print_ligand(*x);
            ++counter;
        }
        
        return counter;
    };

    {
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
          oneapi::tbb::make_filter<mol_vec, std::size_t>(
             oneapi::tbb::filter_mode::serial_out_of_order,
            [&](mol_vec molecules) {
                for (auto& p : molecules) {
                    if (p) {
                        input_queue->enqueue(std::move(p));
                    }
                }
            }
          )
          &
          oneapi::tbb::make_filter<std::size_t, void>(
            oneapi::tbb::filter_mode::parallel,
            [&](std::size_t) {
                drain_ready();
            }
          )
        );

        info("Pipeline done: closing input queue");
        input_queue->close();
        drain_ready();
    }

    // =========================
    // finalize consumer
    // =========================
    // done.store(true);
    // consumer.join();

    info("Output drained to file");
}

} // namespace mudock
    
// namespace mudock {

// // valori “dell’ultima run” nel processo (quindi per-rank)
// static std::atomic<double> g_obs_sum_avg{0.0};
// static std::atomic<std::uint64_t> g_obs_cnt{0};

// void observer_reset() {
//   g_obs_sum_avg.store(0.0, std::memory_order_relaxed);
//   g_obs_cnt.store(0, std::memory_order_relaxed);
// }

// double observer_sum_avg() {
//   return g_obs_sum_avg.load(std::memory_order_relaxed);
// }

// std::uint64_t observer_cnt() {
//   return g_obs_cnt.load(std::memory_order_relaxed);
// }

// } // namespace mudock

// namespace mudock {

//     static inline void print_ligand(const static_molecule& ligand) {
//         std::cout << ligand.properties.get(property_type::NAME) << " "
//                   << ligand.properties.get(property_type::SCORE) << "\n";
//     }

// void run_tbb_pipeline(std::istream& in,
//                       const std::vector<std::string>& configurations,
//                       const knobs& knobs,
//                       genetic_adt_pipeline& pipeline,
//                       std::size_t end,
//                       std::size_t max_tokens,
//                       const double* observer_period_sec,
//                       const double* time_limit_sec) {
//     using mol_vec = parser_filter::mol_vec;

//     auto input_queue  = std::make_shared<safe_stack<static_molecule>>();
//     auto output_queue = std::make_shared<safe_stack<static_molecule>>();

//     std::atomic<std::size_t> processed{0};
//     std::atomic<std::size_t> produced{0};

//     std::atomic<std::size_t> dropped_by_timeout{0};
//     std::atomic<bool> timeout_triggered{false};
//     std::atomic<bool> stop_requested{false};

//     auto drain_ready = [&]() -> std::size_t {
//         std::size_t counter = 0;
//         while (auto x = output_queue->dequeue()) {
//             print_ligand(*x);
//             ++counter;
//         }
//         if (counter) processed.fetch_add(counter, std::memory_order_relaxed);
//         return counter;
//     };

//     mudock::observer_reset();

//     // =========================
//     // Observer thread
//     // =========================
//     std::mutex observer_mutex;
//     std::condition_variable observer_cv;
//     bool observer_stop = false;
//     std::thread observer_thread;

//     const auto start = std::chrono::high_resolution_clock::now();

//     if (observer_period_sec && *observer_period_sec > 0.0) {
//         info("Observer enabled with period: ", *observer_period_sec, " s");
//         observer_thread = std::thread([&]() {
//             std::size_t prev_processed = processed.load(std::memory_order_relaxed);
//             auto prev_time = std::chrono::high_resolution_clock::now();

//             while (true) {
//                 std::unique_lock<std::mutex> lock(observer_mutex);
//                 const bool stop = observer_cv.wait_for(
//                     lock,
//                     std::chrono::duration<double>(*observer_period_sec),
//                     [&]() { return observer_stop; }
//                 );
//                 if (stop) break;
//                 lock.unlock();

//                 // flush output to keep numbers meaningful
//                 drain_ready();

//                 const auto now = std::chrono::high_resolution_clock::now();
//                 const std::size_t now_processed = processed.load(std::memory_order_relaxed);

//                 const std::size_t in_backlog  = input_queue->size();   // requires safe_stack::size()
//                 const std::size_t out_backlog = output_queue->size();  // requires safe_stack::size()

//                 const std::chrono::duration<double> dt = now - prev_time;
//                 const std::size_t delta = now_processed - prev_processed;
//                 const double inst = (dt.count() > 0.0) ? (double(delta) / dt.count()) : 0.0;

//                 const std::chrono::duration<double> total = now - start;
//                 const double avg = (total.count() > 0.0) ? (double(now_processed) / total.count()) : 0.0;

//                 g_obs_sum_avg.fetch_add(avg, std::memory_order_relaxed);
//                 g_obs_cnt.fetch_add(1, std::memory_order_relaxed);

//                 info("Observer: processed=", now_processed,
//                      ", produced=", produced.load(std::memory_order_relaxed),
//                      ", input_backlog=", in_backlog,
//                      ", output_backlog=", out_backlog,
//                      ", inst_throughput=", inst, " ligands/s",
//                      ", avg_throughput=", avg, " ligands/s");

//                 prev_processed = now_processed;
//                 prev_time = now;
//             }
//         });
//     }

//     // =========================
//     // Timer thread (timeout)
//     // =========================
//     std::mutex timer_mutex;
//     std::condition_variable timer_cv;
//     bool timer_cancelled = false;
//     std::thread timer_thread;

//     if (time_limit_sec && *time_limit_sec > 0.0) {
//         info("Time limit enabled: ", *time_limit_sec, " s");
//         timer_thread = std::thread([&]() {
//             std::unique_lock<std::mutex> lock(timer_mutex);
//             const bool cancelled = timer_cv.wait_for(
//                 lock,
//                 std::chrono::duration<double>(*time_limit_sec),
//                 [&]() { return timer_cancelled; }
//             );
//             if (cancelled) return;
//             lock.unlock();

//             // stop producing new chunks
//             stop_requested.store(true, std::memory_order_relaxed);

//             // drop pending work already enqueued (best-effort)
//             dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
//             input_queue->close();

//             timeout_triggered.store(true, std::memory_order_relaxed);

//             info("Time limit reached: discarded ",
//                  dropped_by_timeout.load(std::memory_order_relaxed),
//                  " pending ligand(s) from input queue.");
//         });
//     }

//     {
//         // create workers
//         threadpool pool;
//         manager(configurations, pool, knobs, input_queue, output_queue, pipeline);
//         info("Manager done: workers created");

//         const std::size_t effective_max_tokens = max_tokens == 0 ? 1 : max_tokens;

//         oneapi::tbb::parallel_pipeline(
//             effective_max_tokens,
//             oneapi::tbb::make_filter<void, std::string>(
//                 oneapi::tbb::filter_mode::serial_in_order,
//                 stream_filter(in, end, &stop_requested) // NEW: pass stop
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
//                     std::size_t enq = 0;
//                     for (auto& p : molecules) {
//                         if (!p) continue;

//                         // NEW: stop cooperatively
//                         if (stop_requested.load(std::memory_order_relaxed)) {
//                             // count as dropped logically (optional)
//                             // dropped_by_timeout.fetch_add(1, std::memory_order_relaxed);
//                             continue;
//                         }

//                         input_queue->enqueue(std::move(p));
//                         ++enq;
//                     }
//                     if (enq) produced.fetch_add(enq, std::memory_order_relaxed);
//                     return enq;
//                 }
//             )
//             &
//             oneapi::tbb::make_filter<std::size_t, void>(
//                 oneapi::tbb::filter_mode::serial_out_of_order,
//                 [&](std::size_t) {
//                     drain_ready();
//                 }
//             )
//         );

//         info("Pipeline done: closing input queue");
//         input_queue->close();
//     }

//     // stop observer
//     if (observer_thread.joinable()) {
//         {
//             std::lock_guard<std::mutex> lock(observer_mutex);
//             observer_stop = true;
//         }
//         observer_cv.notify_one();
//         observer_thread.join();
//     }

//     // cancel timer if still waiting
//     if (timer_thread.joinable()) {
//         {
//             std::lock_guard<std::mutex> lock(timer_mutex);
//             timer_cancelled = true;
//         }
//         timer_cv.notify_one();
//         timer_thread.join();
//     }

//     info("Draining output queue");
//     drain_ready();

//     if (timeout_triggered.load(std::memory_order_relaxed)) {
//         info("Dropped ligands due to timeout: ",
//              dropped_by_timeout.load(std::memory_order_relaxed));
//     }
// }

// } // namespace mudock
