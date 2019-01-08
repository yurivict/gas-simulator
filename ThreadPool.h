#pragma once

// based on https://stackoverflow.com/questions/26516683/reusing-thread-in-loop-c

#include <vector>
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

#define THREAD_POOL_DBG(stmt...)

class ThreadPool {
public:
  ThreadPool(int threads) : executing(0), shutdown(false) {
    // Create the specified number of threads
    xthreads.reserve(threads);
    for (int i = 0; i < threads; ++i)
      xthreads.emplace_back(std::bind(&ThreadPool::threadEntry, this, i));
  }
  ~ThreadPool() {
    {
      // Unblock any threads and tell them to stop
      std::unique_lock <std::mutex> l(lock);

      shutdown = true;
      condVarJobAdded.notify_all();
      condVarJobFinished.notify_all();
    }

    // Wait for all threads to stop
    THREAD_POOL_DBG(std::cerr << "Joining threads" << std::endl;)
    for (auto& thread : xthreads)
      thread.join();
  }
  void doJob(std::function <void (void)> func) {
    // Place a job on the queu and unblock a thread
    std::unique_lock<std::mutex> l(lock);

    jobs.emplace(std::move(func));
    condVarJobAdded.notify_one();
  }
  void waitForAll() {
    std::unique_lock<std::mutex> l(lock);

    while (!shutdown && (!jobs.empty() || executing > 0))
      condVarJobFinished.wait(l);
  }
protected:
  void threadEntry(int i) {
    std::function <void (void)> job;

    while (true) {
      {
        std::unique_lock<std::mutex> l(lock);

        while (!shutdown && jobs.empty())
          condVarJobAdded.wait(l);

        if (jobs.empty()) {
          // No jobs to do and we are shutting down
          THREAD_POOL_DBG(std::cerr << "Thread " << i << " terminates" << std::endl;)
          return;
        }

        THREAD_POOL_DBG(std::cerr << "Thread " << i << " does a job" << std::endl;)
        job = std::move(jobs.front());
        jobs.pop();

        ++executing;
      }

      // Do the job without holding any locks
      job();

      { // done
        std::unique_lock<std::mutex> l(lock);
        --executing;
        condVarJobFinished.notify_one();
      }
    }
  }
private:
  std::mutex lock;
  std::condition_variable condVarJobAdded;
  std::condition_variable condVarJobFinished;
  unsigned executing; // how many jobs are currently executing
  bool shutdown;
  std::queue<std::function <void (void)>> jobs;
  std::vector<std::thread> xthreads;
}; // ThreadPool

