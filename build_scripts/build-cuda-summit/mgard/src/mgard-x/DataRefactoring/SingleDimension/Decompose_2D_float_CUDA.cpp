/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/DataRefactoring.hpp"
// clang-format off
namespace mgard_x {

template void decompose_single<2, float, CUDA>(
    Hierarchy<2, float, CUDA> &hierarchy,
    SubArray<2, float, CUDA> &v, SIZE l_target,
    int queue_idx);
} // namespace mgard_x
// clang-format on
