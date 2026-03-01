#pragma once
#include <array>
#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>
#include <atomic>

namespace mudock {

template<class T>
class safe_stack {
public:
  using value_type = T;

private:
  static constexpr std::size_t SHARDS = 16;

  struct shard_t {
    std::vector<std::unique_ptr<value_type>> stack;
    std::mutex m;
  };

  std::array<shard_t, SHARDS> shards_;

  std::condition_variable cv_;
  std::mutex cv_m_;

  std::atomic<bool> closed_{false};
  std::atomic<std::size_t> approx_size_{0};

  static std::size_t tid_hash() {
    return std::hash<std::thread::id>{}(std::this_thread::get_id());
  }

  static std::size_t my_shard() {
    return tid_hash() % SHARDS;
  }

public:
  // non-blocking dequeue: prova tutti gli shard (round-robin)
  [[nodiscard]] inline auto dequeue() {
    std::unique_ptr<value_type> out{};

    const std::size_t start = tid_hash() % SHARDS;
    for (std::size_t k = 0; k < SHARDS; ++k) {
      const std::size_t i = (start + k) % SHARDS;
      auto& sh = shards_[i];
      std::lock_guard<std::mutex> lk(sh.m);
      if (!sh.stack.empty()) {
        out = std::move(sh.stack.back());
        sh.stack.pop_back();
        approx_size_.fetch_sub(1, std::memory_order_relaxed);
        break;
      }
    }
    return out;
  }

  // blocking dequeue: aspetta finché c'è qualcosa o finché è closed
  std::unique_ptr<value_type> dequeue_wait() {
    for (;;) {
      if (auto p = dequeue()) return p;

      // closed + vuoto => fine
      if (closed_.load(std::memory_order_acquire) &&
          approx_size_.load(std::memory_order_relaxed) == 0) {
        return {};
      }

      std::unique_lock<std::mutex> lk(cv_m_);
      cv_.wait(lk, [&]{
        return closed_.load(std::memory_order_acquire) ||
               approx_size_.load(std::memory_order_relaxed) > 0;
      });
    }
  }

  inline std::size_t clear() {
    std::size_t n = 0;

    for (auto& sh : shards_) {
      std::lock_guard<std::mutex> lk(sh.m);
      n += sh.stack.size();
      sh.stack.clear();
    }

    approx_size_.store(0, std::memory_order_relaxed);
    return n;
  }

  void enqueue(std::unique_ptr<value_type> p) {
    if (!p) return;
    if (closed_.load(std::memory_order_acquire)) return;

    const std::size_t i = my_shard();
    {
      auto& sh = shards_[i];
      std::lock_guard<std::mutex> lk(sh.m);
      if (closed_.load(std::memory_order_acquire)) return;
      sh.stack.emplace_back(std::move(p));
    }
    approx_size_.fetch_add(1, std::memory_order_relaxed);
    cv_.notify_one();
  }

  void close() {
    closed_.store(true, std::memory_order_release);
    cv_.notify_all();
  }

  [[nodiscard]] inline std::size_t size() const {
    // approssimato (ma sufficiente per logging/diagnosi)
    return approx_size_.load(std::memory_order_relaxed);
  }
};

} // namespace mudock



// #pragma once

// #include <condition_variable>
// #include <cstdint>
// #include <memory>
// #include <mutex>
// #include <vector>

// namespace mudock {

//   template<class T>
//   class safe_stack {
//   public:
//     using value_type = T;

//   private:
//     std::vector<std::unique_ptr<value_type>> stack;
//     mutable std::mutex mutex;
//     std::condition_variable cv;
//     bool closed = false;

//   public:
//     // non-blocking dequeue, return nullptr if the stack is empty
//     [[nodiscard]] inline auto dequeue() {
//       auto new_element = std::unique_ptr<value_type>{};
//       {
//         std::lock_guard lock{mutex};
//         if (!stack.empty()) {
//           new_element = std::move(stack.back());
//           stack.pop_back();
//         }
//       }
//       return new_element;
//     }

//     // blocking dequeue, return nullptr if the stack is empty and closed
//     std::unique_ptr<value_type> dequeue_wait() {
//       std::unique_lock lock{mutex};
//       // release the lock and reacquire that based on the condition
//       cv.wait(lock, [&]{ return closed || !stack.empty(); });
//       if (stack.empty()) return {};
//       auto p = std::move(stack.back());
//       stack.pop_back();
//       return p;
//     }

//     void enqueue(std::unique_ptr<value_type> p) {
//       if (!p) return;
//       {
//         std::lock_guard lock{mutex};
//         if (closed) return;
//         stack.emplace_back(std::move(p));
//       }
//       cv.notify_one();
//     }

//     // close the stack, i.e., no more element will be enqueued
//     void close() {
//       {
//         std::lock_guard lock{mutex};
//         closed = true;
//       }
//       cv.notify_all();
//     }

//     inline std::size_t clear() {
//       std::lock_guard lock{mutex};
//       const std::size_t n = stack.size();
//       stack.clear();
//       return n;
//     }

//     [[nodiscard]] inline std::size_t size() const {
//       std::lock_guard lock{mutex};
//       return stack.size();
//     }
//   };

// } // namespace mudock
