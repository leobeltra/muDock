#pragma once

#include <istream>
#include <string>
#include <vector>
#include <cstdint>

#include <mudock/mudock.hpp>

namespace mudock {

    // resetta i contatori (da chiamare prima della run)
    void observer_reset();

    // somma degli avg_throughput osservati in questa run (per-rank)
    double observer_sum_avg();

    // numero di osservazioni fatte in questa run (per-rank)
    std::uint64_t observer_cnt();

    static constexpr std::size_t max_tbb_tokens = 4;
    
    // void run_tbb_pipeline(std::istream& in,
    //          const std::vector<std::string>& configurations,
    //          const knobs& knobs,
    //          genetic_adt_pipeline& pipeline, 
    //          std::size_t end = std::numeric_limits<std::size_t>::max(),
    //          std::size_t max_tokens = max_tbb_tokens);
    
    void run_tbb_pipeline(std::istream& in,
                      const std::vector<std::string>& configurations,
                      const knobs& knobs,
                      genetic_adt_pipeline& pipeline,
                      std::size_t end  = std::numeric_limits<std::size_t>::max(),
                      std::size_t max_tokens = max_tbb_tokens,
                      const double* observer_period_sec = nullptr,
                      const double* time_limit_sec = nullptr);
                      
} // namespace mudock
