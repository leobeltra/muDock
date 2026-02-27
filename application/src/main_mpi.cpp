#include <mpi.h>

#include <mudock/mpi_implementation/utilities.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

#include "command_line_args.hpp"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <utility>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  if (rank == 0) {
    mudock::info("[rank 0] nranks = ", nranks);
  }
  mudock::info("[rank ", rank, "] ligand_path = ", args.ligand_path);

  // Pre-check file input (una sola volta)
  {
    std::ifstream probe(args.ligand_path, std::ios::binary);
    if (!probe) {
      mudock::error("[rank ", rank, "] Cannot open ligand file: ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 2);
    }
  }

  // Carico la proteina una volta sola
  auto protein = std::make_shared<mudock::dynamic_molecule>(
      mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  const std::pair<size_t, size_t> range =
      mudock::mpi_splitter_bcast(args.ligand_path, rank, nranks);

  const size_t begin = range.first;
  const size_t end   = range.second;

  // print range per rank
  std::cout << "rank " << rank << ": [" << begin << ", " << end << "]\n";

  // const std::string timing_path = args.timing_path.value_or("mpi_timing.csv");

  // =========================
  // GLOBAL TIMING START
  // // =========================
  // MPI_Barrier(MPI_COMM_WORLD);
  // const double t0 = MPI_Wtime();

    // --- compute our slab ---
  const bool did_work = (begin < end);
  if (did_work) {
    using clock = std::chrono::steady_clock;
    auto start_time  = clock::now();

    std::ifstream in(args.ligand_path, std::ios::binary);
    if (!in) {
      mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);
  // // === splitter ricalcolato ogni iterazione ===
  // const std::pair<size_t, size_t> range =
  //     mudock::mpi_splitter_bcast(args.ligand_path, rank, nranks);

  // const size_t begin = range.first;
  // const size_t end   = range.second;
  // const bool has_work = (begin < end);

  // if (!has_work) {
  //   mudock::info("[rank ", rank, "] No work (empty range).");
  // }

  // if (has_work) {
  //   std::ifstream in(args.ligand_path, std::ios::binary);
  //   if (!in) {
  //     mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
  //     MPI_Abort(MPI_COMM_WORLD, 1);
  //   }

  //   in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);

    mudock::genetic_adt_pipeline pipe{protein};
    mudock::run_tbb_pipeline(in, args.device_confs, args.knobs, pipe, end, 64);

        // per-rank: mean of observer's avg_throughput samples for THIS run
    // const double sum        = mudock::observer_sum_avg();
    // const std::uint64_t cnt = mudock::observer_cnt();
    // local_mean              = (cnt > 0) ? (sum / static_cast<double>(cnt)) : 0.0;

    auto end_time                              = clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    std::cerr << "[GREPME] P" + std::to_string(rank) + "Elapsed time: " + std::to_string(elapsed.count()) +
                     "s \n";

  } else {
    mudock::info("[rank ", rank, "] No work (empty range).");
  }

  // =========================
  // GLOBAL TIMING END
  // =========================
  // MPI_Barrier(MPI_COMM_WORLD);
  // const double t1 = MPI_Wtime();
  // const double local_elapsed = t1 - t0;

  // double max_elapsed = 0.0;
  // double min_elapsed = 0.0;
  // double sum_elapsed = 0.0;

  // MPI_Reduce(&local_elapsed, &max_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  // MPI_Reduce(&local_elapsed, &min_elapsed, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  // MPI_Reduce(&local_elapsed, &sum_elapsed, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  // if (rank == 0) {
  //   const double avg_elapsed = sum_elapsed / static_cast<double>(nranks);

  //   std::ofstream tout(timing_path, std::ios::app);
  //   if (!tout) {
  //     mudock::error("[rank 0] Cannot open timing output file: ", timing_path);
  //   } else {
  //     auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

  //     tout << std::put_time(std::localtime(&now), "%F %T")
  //           << ",nranks=" << nranks
  //           // << ",iter=" << iter
  //           // << ",warmup=" << is_warmup
  //           << ",total_s=" << std::setprecision(10) << max_elapsed
  //           << ",min_s="   << min_elapsed
  //           << ",avg_s="   << avg_elapsed
  //           << ",max_s="   << max_elapsed
  //           << "\n";
  //   }
  // }

  MUDOCK_MARKER_CLOSE;
  // MPI_Finalize();
  return 0;
}



// int main(int argc, char** argv) {
//   MPI_Init(&argc, &argv);

//   int rank = 0, nranks = 1;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//   MPI_Comm_size(MPI_COMM_WORLD, &nranks);

//   const auto args = parse_command_line_arguments(argc, argv);

//   MUDOCK_MARKER_INIT;

//   // --- sanity: input ligands file is readable ---
//   {
//     std::ifstream probe(args.ligand_path, std::ios::binary);
//     if (!probe) {
//       mudock::error("[rank ", rank, "] Cannot open ligand file: ", args.ligand_path);
//       MPI_Abort(MPI_COMM_WORLD, 2);
//     }
//   }

//   // --- each rank reads protein once (outside the loop) ---
//   auto protein = std::make_shared<mudock::dynamic_molecule>(
//       mudock::parser<mudock::dynamic_molecule>(args.protein_path));

//   // --- compute per-rank range once ---
//   const std::pair<size_t, size_t> range = mudock::mpi_splitter(args.ligand_path, rank, nranks);
//   const size_t begin = range.first;
//   const size_t end   = range.second;

//   const bool did_work = (begin < end);
//   if (!did_work) {
//     mudock::info("[rank ", rank, "] No work (empty range).");
//   }

//   // ---- multi-run settings ----
//   constexpr int RUNS = 5;     // total runs
//   constexpr int WARMUP = 1;   // discard first run (run index 0)

//   // rank0 collects kept results (strong scaling cluster throughput per run)
//   std::vector<double> kept;
//   if (rank == 0) kept.reserve(RUNS - WARMUP);

//   for (int r = 0; r < RUNS; ++r) {
//     // all ranks start the run together
//     MPI_Barrier(MPI_COMM_WORLD);

//     double local_mean = 0.0;

//     if (did_work) {
//       // IMPORTANT: re-open input stream every run
//       std::ifstream in(args.ligand_path, std::ios::binary);
//       if (!in) {
//         mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
//         MPI_Abort(MPI_COMM_WORLD, 1);
//       }
//       in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);

//       // IMPORTANT: rebuild pipeline every run (fresh state)
//       mudock::genetic_adt_pipeline pipe{protein};

//       // run pipeline (this should internally call mudock::observer_reset() at start of the run)
//       mudock::run_tbb_pipeline(
//           in,
//           args.device_confs,
//           args.knobs,
//           pipe,
//           end,
//           /*max_tokens=*/4,
//           args.observer ? &*args.observer : nullptr,
//           args.time_limit_sec ? &*args.time_limit_sec : nullptr);

//       // per-rank: mean of observer's avg_throughput samples for THIS run
//       const double sum = mudock::observer_sum_avg();
//       const std::uint64_t cnt = mudock::observer_cnt();
//       local_mean = (cnt > 0) ? (sum / static_cast<double>(cnt)) : 0.0;

//     } else {
//       // no work: contribute 0.0 to strong-scaling sum
//       local_mean = 0.0;
//     }

//     // strong scaling: cluster throughput = sum of per-rank means
//     double cluster_throughput = 0.0;
//     MPI_Reduce(&local_mean, &cluster_throughput, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

//     // discard warmup run(s): keep only r>=WARMUP
//     if (rank == 0 && r >= WARMUP) {
//       kept.push_back(cluster_throughput);
//     }

//     MPI_Barrier(MPI_COMM_WORLD);
//   }

//   // ---- finalize: mean + variance over kept runs ----
//   if (rank == 0) {
//     const int n = static_cast<int>(kept.size()); // should be RUNS-WARMUP (=4)

//     double mean = 0.0;
//     for (double x : kept) mean += x;
//     mean = (n > 0) ? (mean / static_cast<double>(n)) : 0.0;

//     // sample variance (unbiased): divide by (n-1)
//     double var = 0.0;
//     if (n > 1) {
//       for (double x : kept) {
//         const double d = x - mean;
//         var += d * d;
//       }
//       var /= static_cast<double>(n - 1);
//     }

//     const double stddev = (var > 0.0) ? std::sqrt(var) : 0.0;

//     const char* nodes_env = std::getenv("SLURM_JOB_NUM_NODES");
//     const int nodes = nodes_env ? std::atoi(nodes_env) : 0;

//     const auto now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

//     std::filesystem::path out = args.timing_path.value_or("observer_strong_scaling.csv");
//     std::ofstream f(out, std::ios::app);

//     // one CSV line per job
//     f << std::put_time(std::localtime(&now), "%F %T")
//       << ",nodes=" << nodes
//       << ",nranks=" << nranks
//       << ",runs=" << RUNS
//       << ",warmup_discarded=" << WARMUP
//       << ",kept=" << n
//       << ",observer_s=" << args.observer.value_or(0.0)
//       << ",time_limit_s=" << args.time_limit_sec.value_or(0.0)
//       << ",cluster_mean=" << std::setprecision(10) << mean
//       << ",cluster_var=" << std::setprecision(10) << var
//       << ",cluster_stddev=" << std::setprecision(10) << stddev
//       << "\n";
//   }

//   MUDOCK_MARKER_CLOSE;
//   MPI_Finalize();
//   return 0;
// }

