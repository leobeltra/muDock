#pragma once

#include <memory>
#include <string_view>
#include <vector>

#include <mudock/molecule.hpp>
#include <mudock/mudock.hpp>

namespace mudock {
    
    class parser_filter {
        public:
            using mol_vec = std::vector<std::unique_ptr<static_molecule>>;

            mol_vec operator()(std::string_view sv) const;
    };

} // namespace mudock
