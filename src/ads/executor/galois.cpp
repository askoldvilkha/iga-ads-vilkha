// SPDX-FileCopyrightText: 2015 - 2023 Marcin Łoś <marcin.los.91@gmail.com>
// SPDX-License-Identifier: MIT

#include "ads/executor/galois.hpp"

#include "ads/config.hpp"

#include <iostream>

#ifdef ADS_USE_GALOIS

#    include <galois/runtime/Statistics.h>
#    include <galois/substrate/ThreadPool.h>

namespace ads {

galois_executor::galois_executor(int threads) {
    galois::runtime::setStatFile("/dev/null");
    std::cout << "[galois] maxUsableThreads = " 
              << galois::substrate::getThreadPool().getMaxUsableThreads() << "\n";
    thread_count(threads);
}

void galois_executor::thread_count(int threads) {
    galois::setActiveThreads(threads);
    auto actual = galois::setActiveThreads(threads);
    std::cout << "[galois_executor] requested=" << threads << " actual=" << actual << std::endl;

}

}  // namespace ads

#endif  // defined(ADS_USE_GALOIS)
