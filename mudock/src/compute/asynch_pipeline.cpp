#include <fstream>
#include <mudock/compute/asynch_pipeline.hpp>
#include <mudock/format/reader.hpp>
#include <mudock/format/supported_format.hpp>
#include <mudock/splitter.hpp>
#include <stdexcept>
#include <string_view>

namespace mudock {
  void read_and_parse(const std::filesystem::path input_file,
                      std::shared_ptr<mudock::safe_queue<mudock::static_molecule>> input_queue) {
    const auto in_format = mudock::parse_supported_format(input_file);
    constexpr_switch<0, mudock::get_num_supported_format(), 1>(
        [&](const auto format_index) {
          // instantiate the splitter
          const auto format = static_cast<mudock::supported_format>(format_index());
          mudock::splitter<mudock::type_of_format<static_cast<mudock::supported_format>(format_index())>>
              split;

          // open the file with the ligand collection
          auto instream = std::ifstream(input_file);
          if (!instream.good()) {
            throw std::runtime_error("Unable to read the ligand file");
          }

          // declare a buffer of 1MB
          static constexpr auto buffer_size = std::size_t{1048576};
          char buffer[buffer_size];

          // read the stream in chuncks
          while (instream) {
            instream.read(buffer, buffer_size);
            auto ligands_description =
                split(std::string_view{buffer, static_cast<std::size_t>(instream.gcount())});

            // actually parse the ligand and enqueue the data structure
            if constexpr (format == mudock::supported_format::ADTMOL2) {
              for (const auto& description: ligands_description) {
                auto ligand = std::make_unique<mudock::static_molecule>(
                    mudock::parser<mudock::supported_format::ADTMOL2, mudock::static_molecule>(description));
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

          // process the last ligand
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
        },
        in_format);

    // if we reach this point, we need to close the input queue since
    // we have no more ligand to instert in the queue
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
