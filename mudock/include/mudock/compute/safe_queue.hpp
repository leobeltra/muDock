#pragma once

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <mutex>

namespace dpipe {

  template<typename T>
  class finite_queue {
    // this is the actual container of the queue
    std::deque<T> buffer;
    std::size_t max_buffer_size;

    // tto signal when to exit
    bool signal_terminate;

    // these variables avoids busy waiting and actually prevent the thread to consume CPU cycles
    mutable std::mutex queue_mutex;
    std::condition_variable inwork_available;
    std::condition_variable outwork_available;

  public:
    using value_type = T;

    finite_queue(void): max_buffer_size(1), signal_terminate(false) {}

    // this queue cannot be copied or moved around
    finite_queue(const finite_queue &) = delete;
    finite_queue(finite_queue &&)      = delete;

    inline void initialize(const std::size_t max_queue_size) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      max_buffer_size = max_queue_size;
    }

    inline std::size_t size(void) const {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return buffer.size();
    }

    inline std::size_t max_size(void) const { return max_buffer_size; }

    inline bool empty(void) const {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return buffer.empty();
    }

    inline void send_terminate_signal(void) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      signal_terminate = true;
      inwork_available.notify_all();
      outwork_available.notify_all();
    }

    inline void clear_terminate_signal(void) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      signal_terminate = false;
    }

    inline bool get_terminate_signal(void) {
      std::unique_lock<std::mutex> lock(queue_mutex);
      return signal_terminate;
    }

    // this method attempt to remove and return an element from the queue. The output parameters tells if the
    // returned pointer is meaningful or a null pointer.
    // NOTE: this method might block the execution of the application.
    value_type dequeue(bool &is_retrieved) {
      // wait until there is some event
      std::unique_lock<std::mutex> lock(queue_mutex);
      while (!signal_terminate && buffer.empty()) { // spourious events might happens!
        outwork_available.wait(lock);               // release lock -> wait for a wake_up -> reaquire lock
      }

      // if the terminate signal is set, return a nullpointer
      if (signal_terminate) {
        is_retrieved = false;
        return value_type{};
      }

      // otherwise, return the first element available
      auto output_data = std::move(buffer.back());
      buffer.pop_back();
      inwork_available.notify_one();
      is_retrieved = true;
      return output_data;
    }

    // this method attempt to insert an element in the queue. The output flag tells if the input element has
    // been succesfully inserted in the queue. In this case the queue owns the element and the user is not
    // allowed to dereference the pointer. If the operation fails, i.e. the termination signal is set, the owner
    // of the data is still the caller
    void enqueue(value_type &input_data, bool &is_stored) {
      // try to enqueue the element (if there is enough space)
      std::unique_lock<std::mutex> lock(queue_mutex);
      while (!signal_terminate && (buffer.size() >= max_buffer_size)) { // spourious events might happens!
        inwork_available.wait(lock); // release lock -> wait for a wake_up -> reaquire lock
      }

      // if the terminate signal is set, do not enqueue the element
      if (signal_terminate) {
        is_stored = false;
        return;
      }

      // otherwise, store the data in the back of the container
      buffer.emplace_front(std::move(input_data));
      outwork_available.notify_one();
      is_stored = true;
    }
  };
} // namespace dpipe
