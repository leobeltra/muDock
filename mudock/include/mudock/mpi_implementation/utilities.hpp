#pragma once

#include <string>
#include <mudock/format/adt_mol2.hpp>

namespace mudock {

    // helper function to return the offset of the previous marker token
    static std::size_t find_previous_marker(std::istream& file, std::size_t pos) {
        if (pos == 0) return 0;

        const std::string_view marker = adt_mol2_tokens::MOLECULE_TOKEN;
        constexpr std::size_t search_block_size = 2000;

        std::string buffer(search_block_size, '\0');

        while (pos > 0) {
            const std::size_t read_size = std::min(search_block_size, pos);
            const std::size_t start = pos - read_size;

            file.clear();
            file.seekg(static_cast<std::streamoff>(start), std::ios::beg);
            file.read(buffer.data(), static_cast<std::streamsize>(read_size));

            const std::size_t got = static_cast<std::size_t>(file.gcount());
            if (got == 0) break;

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

    // function to split a file into ranges for MPI ranks, aligned to molecule markers
    static std::pair<size_t, size_t>
    mpi_splitter(const std::string& path,
                size_t rank_id,
                size_t num_ranks)
    {
        if (num_ranks == 0 || rank_id >= num_ranks)
            return {0, 0};

        const auto file_size = std::filesystem::file_size(path);
        if (file_size == 0) return {0, 0};
        if (num_ranks == 1) return {0, file_size};

        const size_t chunk_size = file_size / num_ranks;
        size_t start_pos = rank_id * chunk_size;
        size_t end_pos   = (rank_id == num_ranks - 1)
                        ? file_size
                        : start_pos + chunk_size;

        std::ifstream file(path, std::ios::binary);
        if (!file)
            throw std::runtime_error("Failed to open file: " + path);

        // align start to the previous marker (except for rank 0)
        if (rank_id > 0) {
            start_pos = find_previous_marker(file, start_pos);
        }

        // align end to the previous marker (except for the last rank)
        if (rank_id < num_ranks - 1) {
            end_pos = find_previous_marker(file, end_pos);
        }
        
        return {start_pos, end_pos};
    }

} // namespace mudock
