/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Correction/CalcCorrectionND.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrectionND<2, float, CUDA>(
    Hierarchy<2, float, CUDA> &hierarchy,
    SubArray<2, float, CUDA> dcoeff,
    SubArray<2, float, CUDA> &dcorrection, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
