/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Correction/CalcCorrection3D.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrection3D<1, float, SERIAL>(
    Hierarchy<1, float, SERIAL> &hierarchy,
    SubArray<1, float, SERIAL> dcoeff,
    SubArray<1, float, SERIAL> &dcorrection, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
