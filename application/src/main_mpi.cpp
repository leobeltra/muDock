#include "command_line_args.hpp"

#include <mpi.h>
#include <mudock/mpi_implementation/utilities.hpp>
#include <mudock/tbb_implementation/tbb_pipeline.hpp>
#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  int rank = 0, nranks = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);

  const auto args = parse_command_line_arguments(argc, argv);

  MUDOCK_MARKER_INIT;

  // each rank reads the protein
  mudock::info("[rank ", rank, "] ligand_path = ", args.ligand_path);

  std::ifstream probe(args.ligand_path, std::ios::binary);
  if (!probe) {
    mudock::error("[rank ", rank, "] Cannot open ligand file: ", args.ligand_path);
    MPI_Abort(MPI_COMM_WORLD, 2);
  }
  
  // optional: print size via tellg to confirm
  // probe.seekg(0, std::ios::end);
  // mudock::info("[rank ", rank, "] ligand file size = ", (std::size_t)probe.tellg());
  // probe.close();  
  auto protein =
        std::make_shared<mudock::dynamic_molecule>(
            mudock::parser<mudock::dynamic_molecule>(args.protein_path));

  // MPI splitter returns the range as a std::pair
  std::pair<size_t, size_t> range =
      mudock::mpi_splitter(args.ligand_path, rank, nranks);

  const size_t begin = range.first;
  const size_t end   = range.second;

  // If a rank does not find any work, return
  if (begin >= end) {
    mudock::info("[rank ", rank, "] No work (empty range).");
    MUDOCK_MARKER_CLOSE;
    MPI_Finalize();
    return 0;
  }

  // Debug print to ensure ranges are correct
  // mudock::info("[rank ", rank, "] Ligand range: [", begin, ", ", end, ")");

  std::ifstream in(args.ligand_path, std::ios::binary);
  if (!in) {
    mudock::error("[rank ", rank, "] Can't open input file ", args.ligand_path);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  // each rank seeks to its begin
  in.seekg(static_cast<std::streamoff>(begin), std::ios::beg);

  mudock::genetic_adt_pipeline pipe{protein};

  // stream limited to end
  mudock::run_tbb_pipeline(in, args.device_confs, args.knobs, pipe, end);

  MUDOCK_MARKER_CLOSE;
  MPI_Finalize();
  return 0;
}
