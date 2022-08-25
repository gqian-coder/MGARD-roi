/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/CopyND/CopyND.hpp"
// clang-format off
namespace mgard_x {

template void CopyND<5, double, CUDA>(
    SubArray<5, double, CUDA> dinput,
    SubArray<5, double, CUDA> &doutput, int queue_idx);

} // namespace mgard_x
// clang-format on
