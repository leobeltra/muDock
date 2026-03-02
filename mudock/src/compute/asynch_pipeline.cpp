#include <fstream>
#include <mpi.h>
#include <mudock/compute/asynch_pipeline.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/mpi_implementation/utilities.hpp>
#include <mudock/splitter.hpp>
#include <stdexcept>
#include <string>
#include <string_view>

static std::pair<size_t, size_t> slab_finder(size_t file_size, size_t rank, size_t comm_size) {
  if (comm_size == 0 || rank >= comm_size)
    return {0, 0};
  if (file_size == 0)
    return {0, 0};
  if (comm_size == 1)
    return {0, file_size};

  const size_t chunk_size = file_size / comm_size;
  const auto begin        = rank * chunk_size;
  return {begin, (rank == comm_size - 1) ? file_size : begin + chunk_size};
}

namespace mudock {
  void read_and_parse(const std::filesystem::path input_file,
                      std::shared_ptr<mudock::safe_queue<mudock::static_molecule>> input_queue) {
    // figure out our place in the universe
    int rank = 0, size = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // get the file size
    unsigned long long file_size = 0;
    if (rank == 0) {
      file_size = std::filesystem::file_size(input_file);
    }
    MPI_Bcast(&file_size, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

    // find the slab that we need to compute
    auto [begin, end] = slab_finder(file_size, rank, size);

    // open the file with the ligand collection
    auto instream = std::ifstream(input_file);
    if (!instream.good()) {
      throw std::runtime_error("Unable to read the ligand file");
    }

    // align the begin and end of the slab according to the content
    begin = align_marker(instream, begin);

    // print debug message
    std::cerr << "[GREPME] P" + std::to_string(rank) + " - Slab size: " + std::to_string(end - begin) + " \n";

    // make sure that we need to compute something
    if (begin < end) {
      // set the file stream to the initial offset for our slab
      instream.seekg(static_cast<std::streamoff>(begin), std::ios::beg);
      auto current_slab_offset = begin;

      // start to read and parse the file
      const auto in_format = mudock::parse_supported_format(input_file);
      constexpr_switch<0, mudock::get_num_supported_format(), 1>(
          [&](const auto format_index) {
            // instantiate the splitter
            const auto format = static_cast<mudock::supported_format>(format_index());
            mudock::splitter<mudock::type_of_format<static_cast<mudock::supported_format>(format_index())>>
                split;

            // declare a buffer of 1MB
            static constexpr auto buffer_size = std::size_t{1048576};
            char buffer[buffer_size];

            // read the stream in chuncks
            mudock::info("Starting to read and parse ligands ...");
            auto need_to_read = true;
            while (need_to_read && instream) {
              // get all the descriptions in the read batch
              instream.read(buffer, buffer_size);
              auto ligands_description =
                  split(std::string_view{buffer, static_cast<std::size_t>(instream.gcount())});

              // make sure to compute only the molecules assigned to us
              auto valid_size = std::size_t{0};
              for (const auto& desc: ligands_description) {
                if (current_slab_offset < end) {
                  valid_size += 1;
                  current_slab_offset += desc.size();
                } else {
                  need_to_read = false;
                  break;
                }
              }
              ligands_description.resize(valid_size);

              // actually parse the ligand and enqueue the data structure
              if constexpr (format == mudock::supported_format::ADTMOL2) {
                for (const auto& description: ligands_description) {
                  auto ligand = std::make_unique<mudock::static_molecule>(
                      mudock::parser<mudock::supported_format::ADTMOL2, mudock::static_molecule>(
                          description));
                  auto is_stored = false;
                  input_queue->enqueue(ligand, is_stored);
                  assert(is_stored);
                }
              } else {
                for (const auto& description: ligands_description) {
                  try {
                    auto ligand = std::make_unique<mudock::static_molecule>(
                        mudock::parser<format, mudock::static_molecule>(description));
                    auto is_stored = false;
                    input_queue->enqueue(ligand, is_stored);
                    assert(is_stored);
                  } catch (...) {}
                }
              }
            }

            // process the last ligand (if we need to process it)
            if (need_to_read) {
              const auto description = split.flush();
              if constexpr (format == mudock::supported_format::ADTMOL2) {
                auto ligand = std::make_unique<mudock::static_molecule>(
                    mudock::parser<mudock::supported_format::ADTMOL2, mudock::static_molecule>(description));
                auto is_stored = false;
                input_queue->enqueue(ligand, is_stored);
                assert(is_stored);

              } else {
                try {
                  auto ligand = std::make_unique<mudock::static_molecule>(
                      mudock::parser<format, mudock::static_molecule>(description));
                  auto is_stored = false;
                  input_queue->enqueue(ligand, is_stored);
                  assert(is_stored);
                } catch (...) {}
              }
            }
          },
          in_format);
    }

    // if we reach this point, we need to close the input queue since
    // we have no more ligand to instert in the queue
    mudock::info("All the ligands has been read and parsed!");
    input_queue->send_terminate_signal();
  }

  void write(std::shared_ptr<mudock::safe_queue<mudock::static_molecule>> output_queue) {
    bool is_retrieved = false;
    for (auto ligand = output_queue->dequeue(is_retrieved); ligand;
         ligand      = output_queue->dequeue(is_retrieved)) {
      if (is_retrieved)
        std::cout << ligand->properties.get(mudock::property_type::NAME) + " " +
                         ligand->properties.get(mudock::property_type::SCORE) + "\n";
    }
  }

} // namespace mudock
