/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/DataRefactoring.hpp"
// clang-format off
namespace mgard_x {

template void decompose_single<5, float, SERIAL>(
    Hierarchy<5, float, SERIAL> &hierarchy,
    SubArray<5, float, SERIAL> &v, SIZE l_target,
    int queue_idx);
} // namespace mgard_x
// clang-format on
