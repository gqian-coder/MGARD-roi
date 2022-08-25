/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/CopyND/CopyND.hpp"
// clang-format off
namespace mgard_x {

template void CopyND<3, float, CUDA>(
    SubArray<3, float, CUDA> dinput,
    SubArray<3, float, CUDA> &doutput, int queue_idx);

} // namespace mgard_x
// clang-format on
