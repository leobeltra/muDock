#pragma once

#include <filesystem>
#include <mudock/format/adt_mol2.hpp>
#include <string>

namespace mudock {

  // helper function to return the offset of the previous marker token
  static std::size_t align_marker(std::istream& file, std::size_t pos) {
    const std::string_view marker = adt_mol2_tokens::MOLECULE_TOKEN;

    // define the reading buffer that we go back to search for the starting molecule
    constexpr std::size_t search_block_size = 1024 * 40; // 40Kb
    std::string buffer(search_block_size, '\0');

    // read the file to look for the correct molecule
    file.seekg(static_cast<std::streamoff>(pos), std::ios::beg);
    file.read(buffer.data(), static_cast<std::streamsize>(search_block_size));
    buffer.resize(file.gcount());

    // find the first occurrence of the maker (if any)
    const auto first_occur = buffer.find(marker);
    return first_occur != std::string::npos ? pos + first_occur : std::string::npos;
  }

} // namespace mudock
