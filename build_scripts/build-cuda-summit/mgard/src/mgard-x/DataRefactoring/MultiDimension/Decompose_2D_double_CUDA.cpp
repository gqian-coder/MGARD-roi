/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/DataRefactoring.hpp"
// clang-format off
namespace mgard_x {

template void decompose<2, double, CUDA>(
    Hierarchy<2, double, CUDA> &hierarchy,
    SubArray<2, double, CUDA> &v, SIZE l_target,
    int queue_idx);
} // namespace mgard_x
// clang-format on
