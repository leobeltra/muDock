#include "command_line_args.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <mpi.h>
#include <mudock/molecule.hpp>
#include <mudock/mpi_implementation/utilities.hpp>
#include <mudock/mudock.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <vector>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  // --- each rank reads protein once (outside the loop) ---
  auto protein =
      std::make_shared<mudock::dynamic_molecule>(mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  // --- compute per-rank range once ---
  const std::pair<size_t, size_t> range = mudock::mpi_splitter(args.ligand_path, rank, nranks);
  const size_t begin                    = range.first;
  const size_t end                      = range.second;

  // --- compute our slab ---
  const bool did_work = (begin < end);
  if (did_work) {
    using clock = std::chrono::steady_clock;
    auto start  = clock::now();

    std::ifstream in(args.ligand_path, std::ios::binary);
    if (!in) {
      mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);

    // create the pipeline
    mudock::genetic_adt_pipeline pipe{protein};

    // run pipeline (this should internally call mudock::observer_reset() at start of the run)
    mudock::run_tbb_pipeline(in,
                             args.device_confs,
                             args.knobs,
                             pipe,
                             end,
                             /*max_tokens=*/4);

    auto stop                             = clock::now();
    std::chrono::duration<double> elapsed = stop - start;

    std::cerr << "[GREPME] P" + std::to_string(rank) + "Elapsed time: " + std::to_string(elapsed.count()) +
                     "s \n";

  } else {
    mudock::info("[rank ", rank, "] No work (empty range).");
  }

  MUDOCK_MARKER_CLOSE;
  MPI_Finalize();
  return 0;
}
