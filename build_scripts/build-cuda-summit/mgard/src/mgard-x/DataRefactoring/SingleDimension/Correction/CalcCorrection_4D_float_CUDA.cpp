/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/Correction/CalcCorrection.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrection<4, float, CUDA>(
    Hierarchy<4, float, CUDA> &hierarchy,
    SubArray<4, float, CUDA> &coeff,
    SubArray<4, float, CUDA> &correction,
    SIZE curr_dim, SIZE l, int queue_idx);

} // namespace mgard_x
// clang-format off
