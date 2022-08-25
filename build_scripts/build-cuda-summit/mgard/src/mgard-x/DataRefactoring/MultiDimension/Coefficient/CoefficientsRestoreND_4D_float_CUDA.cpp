/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Coefficient/CoefficientsRestoreND.hpp"
// clang-format off
namespace mgard_x {

template void CoefficientsRestoreND<4, float, CUDA>(
    Hierarchy<4, float, CUDA> &hierarchy,
    SubArray<4, float, CUDA> dinput1,
    SubArray<4, float, CUDA> dinput2,
    SubArray<4, float, CUDA> &doutput, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
