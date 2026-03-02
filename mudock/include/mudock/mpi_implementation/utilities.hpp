#pragma once

#include <filesystem>
#include <mudock/format/adt_mol2.hpp>
#include <string>

namespace mudock {

  // helper function to return the offset of the previous marker token
  static std::size_t find_previous_marker(std::istream& file, std::size_t pos) {
    if (pos == 0)
      return 0;

    const std::string_view marker           = adt_mol2_tokens::MOLECULE_TOKEN;
    constexpr std::size_t search_block_size = 2000;

    std::string buffer(search_block_size, '\0');

    while (pos > 0) {
      const std::size_t read_size = std::min(search_block_size, pos);
      const std::size_t start     = pos - read_size;

      file.clear();
      file.seekg(static_cast<std::streamoff>(start), std::ios::beg);
      file.read(buffer.data(), static_cast<std::streamsize>(read_size));

      const std::size_t got = static_cast<std::size_t>(file.gcount());
      if (got == 0)
        break;

      const std::string_view view(buffer.data(), got);
      const std::size_t cut = view.rfind(marker);

      // if marker is found, return its offset
      if (cut != std::string_view::npos) {
        return start + cut;
      }

      // go to the previous block
      pos = start;
    }

    return 0;
  }

} // namespace mudock
