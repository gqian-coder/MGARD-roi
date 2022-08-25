/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/Correction/CalcCorrection.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrection<5, double, SERIAL>(
    Hierarchy<5, double, SERIAL> &hierarchy,
    SubArray<5, double, SERIAL> &coeff,
    SubArray<5, double, SERIAL> &correction,
    SIZE curr_dim, SIZE l, int queue_idx);

} // namespace mgard_x
// clang-format off
