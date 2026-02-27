#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <memory>
#include <mudock/compute/manager.hpp>
#include <mudock/compute/safe_stack.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>
#include <mudock/tbb_implementation/parser_filter.hpp>
#include <mudock/tbb_implementation/stream_filter.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mutex>
#include <oneapi/tbb/parallel_pipeline.h>
#include <thread>
#include <vector>

namespace mudock {

  // valori “dell’ultima run” nel processo (quindi per-rank)
  static std::atomic<double> g_obs_sum_avg{0.0};
  static std::atomic<std::uint64_t> g_obs_cnt{0};

  void observer_reset() {
    g_obs_sum_avg.store(0.0, std::memory_order_relaxed);
    g_obs_cnt.store(0, std::memory_order_relaxed);
  }

  double observer_sum_avg() { return g_obs_sum_avg.load(std::memory_order_relaxed); }

  std::uint64_t observer_cnt() { return g_obs_cnt.load(std::memory_order_relaxed); }

} // namespace mudock

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
                        std::size_t max_tokens,
                        const double* observer_period_sec,
                        const double* time_limit_sec) {
    using mol_vec = parser_filter::mol_vec;

    auto input_queue  = std::make_shared<safe_stack<static_molecule>>();
    auto output_queue = std::make_shared<safe_stack<static_molecule>>();

    std::atomic<std::size_t> processed{0};
    std::atomic<std::size_t> produced{0};

    std::atomic<bool> timeout_triggered{false};
    std::atomic<bool> stop_requested{false};

    auto drain_ready = [&]() -> std::size_t {
      std::size_t counter = 0;
      while (auto x = output_queue->dequeue()) {
        print_ligand(*x);
        ++counter;
      }
      if (counter)
        processed.fetch_add(counter, std::memory_order_relaxed);
      return counter;
    };

    {
      // create workers
      threadpool pool;
      manager(configurations, pool, knobs, input_queue, output_queue, pipeline);
      info("Manager done: workers created");

      const std::size_t effective_max_tokens = max_tokens == 0 ? 1 : max_tokens;

      oneapi::tbb::parallel_pipeline(
          effective_max_tokens,
          oneapi::tbb::make_filter<void, std::string>(
              oneapi::tbb::filter_mode::serial_in_order,
              stream_filter(in, end, &stop_requested) // NEW: pass stop
              ) &
              oneapi::tbb::make_filter<std::string, mol_vec>(oneapi::tbb::filter_mode::parallel,
                                                             parser_filter()) &
              oneapi::tbb::make_filter<mol_vec, std::size_t>(oneapi::tbb::filter_mode::serial_out_of_order,
                                                             [&](mol_vec molecules) -> std::size_t {
                                                               for (auto& p: molecules) {
                                                                 if (!p)
                                                                   continue;

                                                                 // NEW: stop cooperatively
                                                                 if (stop_requested.load(
                                                                         std::memory_order_relaxed)) {
                                                                   // count as dropped logically (optional)
                                                                   // dropped_by_timeout.fetch_add(1, std::memory_order_relaxed);
                                                                   continue;
                                                                 }

                                                                 input_queue->enqueue(std::move(p));
                                                               }
                                                             }) &
              oneapi::tbb::make_filter<std::size_t, void>(oneapi::tbb::filter_mode::serial_out_of_order,
                                                          [&]() { drain_ready(); }));

      info("Pipeline done: closing input queue");
      input_queue->close();
    }

    info("Draining output queue");
    drain_ready();
  }

} // namespace mudock
