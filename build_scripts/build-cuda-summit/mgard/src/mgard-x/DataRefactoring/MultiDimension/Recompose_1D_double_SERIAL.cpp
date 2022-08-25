/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/DataRefactoring.hpp"
// clang-format off
namespace mgard_x {

template void recompose<1, double, SERIAL>(
    Hierarchy<1, double, SERIAL> &hierarchy,
    SubArray<1, double, SERIAL> &v, SIZE l_target,
    int queue_idx);
} // namespace mgard_x
// clang-format on
