#pragma once

#include <istream>
#include <string>
#include <vector>
#include <cstdint>

#include <mudock/mudock.hpp>

namespace mudock {

    static constexpr std::size_t max_tbb_tokens = 4;
    
    void run_tbb_pipeline(std::istream& in,
             const std::vector<std::string>& configurations,
             const knobs& knobs,
             genetic_adt_pipeline& pipeline, 
             std::size_t end = std::numeric_limits<std::size_t>::max(),
             std::size_t max_tokens = max_tbb_tokens);
    
    // Dedicated observer thread (if observer_period_sec is set) that periodically prints throughput and backlog info         
    // void observer_reset();

    // double observer_sum_avg();

    // std::uint64_t observer_cnt();

    // void run_tbb_pipeline(std::istream& in,
    //                   const std::vector<std::string>& configurations,
    //                   const knobs& knobs,
    //                   genetic_adt_pipeline& pipeline,
    //                   std::size_t end  = std::numeric_limits<std::size_t>::max(),
    //                   std::size_t max_tokens = max_tbb_tokens,
    //                   const double* observer_period_sec = nullptr,
    //                   const double* time_limit_sec = nullptr);
                      
} // namespace mudock
