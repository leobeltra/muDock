#pragma once

#include <filesystem>
#include <fstream>
#include <memory>
#include <mudock/compute/safe_queue.hpp>
#include <mudock/molecule.hpp>

namespace mudock {

  void read_and_parse(const std::filesystem::path input_file,
                      std::shared_ptr<mudock::safe_queue<mudock::static_molecule>> input_stack);

  void write(std::shared_ptr<mudock::safe_queue<mudock::static_molecule>> output_stack);

} // namespace mudock
