// #pragma once

// #include <istream>
// #include <string>
// #include <oneapi/tbb/parallel_pipeline.h>

// namespace mudock {

//   static constexpr std::size_t max_bytes_per_slice = 20000;

//   class stream_filter {
//     std::istream& stream_;
//     size_t end_;

//   public:
//     explicit stream_filter(std::istream& in, std::size_t end = std::numeric_limits<std::size_t>::max());
    
//     // Source filter: reads from the input stream and produces strings
//     std::string operator()(oneapi::tbb::flow_control& fc, 
//                            std::size_t max_bytes = max_bytes_per_slice) const;
//   };

// } // namespace mudock 

// stream_filter.hpp
#pragma once
#include <oneapi/tbb/parallel_pipeline.h>
#include <atomic>
#include <istream>
#include <string>

namespace mudock {

static constexpr std::size_t max_bytes_per_slice = 20000;

class stream_filter {
  std::istream& stream_;
  std::size_t end_;
  const std::atomic<bool>* stop_{nullptr};

public:
  explicit stream_filter(std::istream& in,
                         std::size_t end = std::numeric_limits<std::size_t>::max(),
                         const std::atomic<bool>* stop = nullptr);

  // questo è quello che TBB chiama
  std::string operator()(oneapi::tbb::flow_control& fc) const {
    return (*this)(fc, max_bytes_per_slice);
  }

  // stesso di prima, ma SENZA default -> niente ambiguità
  std::string operator()(oneapi::tbb::flow_control& fc,
                         std::size_t max_bytes) const;
};

} // namespace mudock
