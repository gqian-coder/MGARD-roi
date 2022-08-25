/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Correction/CalcCorrectionND.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrectionND<1, float, CUDA>(
    Hierarchy<1, float, CUDA> &hierarchy,
    SubArray<1, float, CUDA> dcoeff,
    SubArray<1, float, CUDA> &dcorrection, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
