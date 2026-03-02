#include "command_line_args.hpp"

#include <mpi.h>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <memory>
#include <thread>
#include <mudock/chem/autodock_grid_types.hpp>
#include <mudock/compute/manager.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/likwid_utils.hpp>
#include <mudock/mpi_implementation/utilities.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  // read and parse the target protein
  mudock::info("[rank ", rank, "] Reading and parsing protein ", args.protein_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  // read all the ligands description and split MPI ranges as in main_mpi
  mudock::info("[rank ", rank, "] Reading ligand ", args.ligand_path, " ...");
  const auto in_format = mudock::parse_supported_format(args.ligand_path);
  auto input_queue     = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();

  const auto range =
      mudock::mpi_splitter_bcast(args.ligand_path.string(), rank, nranks, MPI_COMM_WORLD);
  const std::size_t begin = static_cast<std::size_t>(range.first);
  const std::size_t end   = static_cast<std::size_t>(range.second);
  mudock::info("[rank ", rank, "] range: [", begin, ", ", end, "]");

  std::string input_text;
  if (begin < end) {
    std::ifstream in(args.ligand_path, std::ios::binary);
    if (!in) {
      mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);
    const std::size_t bytes = end - begin;
    input_text.resize(bytes);
    in.read(input_text.data(), static_cast<std::streamsize>(bytes));
    input_text.resize(static_cast<std::size_t>(in.gcount()));
  }

  constexpr_switch<0, mudock::get_num_supported_format(), 1>(
      [&](const auto format_index) {
        const auto format = static_cast<mudock::supported_format>(format_index());
        mudock::splitter<mudock::type_of_format<static_cast<mudock::supported_format>(format_index())>> split;
        auto ligands_description = split(std::move(input_text));
        ligands_description.emplace_back(split.flush());

        // parse the input ligands and put them in a stack that we can compute
        mudock::info("[rank ", rank, "] Parsing ", ligands_description.size(), " ligand(s) ...");
        if constexpr (format == mudock::supported_format::ADTMOL2) {
#ifdef _OPENMP
#pragma omp parallel for shared(input_queue)
#endif
          for (const auto& description: ligands_description) {
            auto ligand = std::make_unique<mudock::static_molecule>(
                mudock::parser<mudock::supported_format::ADTMOL2, mudock::static_molecule>(description));
            input_queue->enqueue(std::move(ligand));
          }
        } else {
          for (const auto& description: ligands_description) {
            try {
              auto ligand = std::make_unique<mudock::static_molecule>(
                  mudock::parser<format, mudock::static_molecule>(description));
              input_queue->enqueue(std::move(ligand));
            } catch (...) {}
          }
        }
      },
      in_format);

  // compute all the ligands according to the input configuration
  mudock::info("[rank ", rank, "] Virtual screening the ligands ...");
  mudock::genetic_adt_pipeline pipe{protein};

  auto output_queue = std::make_shared<mudock::safe_stack<mudock::static_molecule>>();
  const auto start  = std::chrono::high_resolution_clock::now();
  std::atomic<std::size_t> dropped_by_timeout{0};
  std::atomic<bool> timeout_triggered{false};
  {
    std::mutex observer_mutex;
    std::condition_variable observer_cv;
    bool observer_stop = false;
    std::thread observer_thread;
    if (args.observer && *args.observer > 0.0) {
      mudock::info("[rank ", rank, "] Observer enabled with period: ", *args.observer, " s");
      observer_thread = std::thread([&]() {
        std::size_t prev_processed = output_queue->size();
        auto prev_time             = std::chrono::high_resolution_clock::now();
        while (true) {
          std::unique_lock<std::mutex> lock(observer_mutex);
          const bool stop = observer_cv.wait_for(lock,
                                                 std::chrono::duration<double>(*args.observer),
                                                 [&]() { return observer_stop; });
          if (stop) {
            break;
          }
          lock.unlock();

          const auto now                  = std::chrono::high_resolution_clock::now();
          const std::size_t now_processed = output_queue->size();
          const std::size_t in_backlog    = input_queue->size();

          const std::chrono::duration<double> dt = now - prev_time;
          const std::size_t delta_processed       = now_processed - prev_processed;
          const double inst_throughput =
              dt.count() > 0.0 ? static_cast<double>(delta_processed) / dt.count() : 0.0;
          const std::chrono::duration<double> total = now - start;
          const double avg_throughput =
              total.count() > 0.0 ? static_cast<double>(now_processed) / total.count() : 0.0;

          mudock::info("[rank ",
                       rank,
                       "] Observer: processed=",
                       now_processed,
                       ", input_backlog=",
                       in_backlog,
                       ", inst_throughput=",
                       inst_throughput,
                       " ligands/s, avg_throughput=",
                       avg_throughput,
                       " ligands/s");

          prev_processed = now_processed;
          prev_time      = now;
        }
      });
    }

    std::mutex timer_mutex;
    std::condition_variable timer_cv;
    bool timer_cancelled = false;
    std::thread timer_thread;
    if (args.time_limit_sec && *args.time_limit_sec > 0.0) {
      mudock::info("[rank ", rank, "] Time limit enabled: ", *args.time_limit_sec, " s");
      timer_thread = std::thread([&]() {
        std::unique_lock<std::mutex> lock(timer_mutex);
        const bool cancelled = timer_cv.wait_for(lock,
                                                 std::chrono::duration<double>(*args.time_limit_sec),
                                                 [&]() { return timer_cancelled; });
        if (cancelled) {
          return;
        }
        lock.unlock();
        dropped_by_timeout.store(input_queue->clear(), std::memory_order_relaxed);
        timeout_triggered.store(true, std::memory_order_relaxed);
        mudock::info("[rank ",
                     rank,
                     "] Time limit reached: discarded ",
                     dropped_by_timeout.load(std::memory_order_relaxed),
                     " pending ligand(s) from input queue.");
      });
    }

    {
      auto threadpool = mudock::threadpool();
      mudock::manager(args.device_confs, threadpool, args.knobs, input_queue, output_queue, pipe);
      input_queue->close(); // signal that no more ligand will be enqueued, so the workers can stop when they finish the backlog
      mudock::info("[rank ", rank, "] All workers have been created!");
    } // threadpool destructor waits for workers; computation is complete here

    mudock::info("Computation complete, shutting down observer and timer threads ...");

    if (observer_thread.joinable()) {
      {
        std::lock_guard<std::mutex> lock(observer_mutex);
        observer_stop = true;
      }
      observer_cv.notify_one();
      observer_thread.join();
    }

    if (timer_thread.joinable()) {
      {
        std::lock_guard<std::mutex> lock(timer_mutex);
        timer_cancelled = true;
      }
      timer_cv.notify_one();
      timer_thread.join();
    }
  } // when we exit from this block the computation is complete

  if (timeout_triggered.load(std::memory_order_relaxed)) {
    mudock::info("[rank ",
                 rank,
                 "] Dropped ligands due to timeout: ",
                 dropped_by_timeout.load(std::memory_order_relaxed));
  }

  // after the computation it will be nice to print the score of all the molecules
  mudock::info("[rank ", rank, "] Printing the scores ...");
  for (auto ligand = output_queue->dequeue(); ligand; ligand = output_queue->dequeue()) {
    std::cout << ligand->properties.get(mudock::property_type::NAME) << " "
              << ligand->properties.get(mudock::property_type::SCORE) << std::endl;
  }

  MUDOCK_MARKER_CLOSE;
  MPI_Finalize();

  // if we reach this statement we completed successfully the run
  mudock::info("[rank ", rank, "] All Done!");
  return EXIT_SUCCESS;
}
