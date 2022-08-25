/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/Correction/CalcCorrection.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrection<3, float, CUDA>(
    Hierarchy<3, float, CUDA> &hierarchy,
    SubArray<3, float, CUDA> &coeff,
    SubArray<3, float, CUDA> &correction,
    SIZE curr_dim, SIZE l, int queue_idx);

} // namespace mgard_x
// clang-format off
