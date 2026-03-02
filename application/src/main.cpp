#include "command_line_args.hpp"

#include <atomic>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <mudock/mudock.hpp>
#include <mutex>
#include <thread>

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  // start measuring the walltime
  using clock = std::chrono::high_resolution_clock;
  auto start  = clock::now();

  // read and parse the target protein
  mudock::info("Reading and parsing protein ", args.protein_path, " ...");
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  // declare and initialize the data structures that are shared among workers
  auto input_queue  = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  auto output_queue = std::make_shared<mudock::safe_queue<mudock::static_molecule>>();
  input_queue->initialize(1000);
  output_queue->initialize(1000);
  mudock::genetic_adt_pipeline pipe{protein};

  // start the reader
  mudock::info("Reading ligand in \"", args.ligand_path, "\" ...");
  auto reader = std::async(std::launch::async, mudock::read_and_parse, args.ligand_path, input_queue);

  // start the writer
  auto writer = std::async(std::launch::async, mudock::write, output_queue);

  // start the observer
  std::mutex observer_mutex;
  std::condition_variable observer_cv;
  bool observer_stop = false;
  std::thread observer_thread;
  if (args.observer && *args.observer > 0.0) {
    mudock::info("Observer enabled with period: ", *args.observer, " s");
    observer_thread = std::thread([&]() {
      std::size_t prev_processed = output_queue->size();
      auto prev_time             = std::chrono::high_resolution_clock::now();
      while (true) {
        std::unique_lock<std::mutex> lock(observer_mutex);
        const bool stop = observer_cv.wait_for(lock, std::chrono::duration<double>(*args.observer), [&]() {
          return observer_stop;
        });
        if (stop) {
          break;
        }
        lock.unlock();

        const auto now                  = std::chrono::high_resolution_clock::now();
        const std::size_t now_processed = output_queue->get_global_counter();
        const std::size_t in_backlog    = input_queue->get_global_counter();

        const std::chrono::duration<double> dt = now - prev_time;
        const std::size_t delta_processed      = now_processed - prev_processed;
        const double inst_throughput =
            dt.count() > 0.0 ? static_cast<double>(delta_processed) / dt.count() : 0.0;
        const std::chrono::duration<double> total = now - start;
        const double avg_throughput =
            total.count() > 0.0 ? static_cast<double>(now_processed) / total.count() : 0.0;

        mudock::info("Observer: processed=",
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

  // spaw
  mudock::info("Spawing the working threads ...");
  auto threadpool = mudock::threadpool();
  mudock::manager(args.device_confs, threadpool, args.knobs, input_queue, output_queue, pipe);

  // wait until the reader and the pipeline completed the execution
  reader.wait();
  threadpool.wait();

  // send the terminate message to the output queue
  output_queue->send_terminate_signal();

  // wait until the writer complete its job
  writer.wait();
  mudock::info("The computation is done!");

  // write a message with the measured walltime
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  auto stop                             = clock::now();
  std::chrono::duration<double> elapsed = stop - start;
  std::cerr << "[GREPME] P" + std::to_string(rank) + " - Elapsed time: " + std::to_string(elapsed.count()) +
                   "s (processed " + std::to_string(output_queue->get_global_counter()) + " ligands) \n";
  // stop the obseerver
  if (observer_thread.joinable()) {
    {
      std::lock_guard<std::mutex> lock(observer_mutex);
      observer_stop = true;
    }
    observer_cv.notify_one();
    observer_thread.join();
  }

  MUDOCK_MARKER_CLOSE;

  // if we reach this statement we completed successfully the run
  mudock::info("All Done!");
  MPI_Finalize();
  return EXIT_SUCCESS;
}
