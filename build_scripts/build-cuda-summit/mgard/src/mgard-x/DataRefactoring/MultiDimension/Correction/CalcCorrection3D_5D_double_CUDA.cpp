/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Correction/CalcCorrection3D.hpp"
// clang-format off
namespace mgard_x {

template void CalcCorrection3D<5, double, CUDA>(
    Hierarchy<5, double, CUDA> &hierarchy,
    SubArray<5, double, CUDA> dcoeff,
    SubArray<5, double, CUDA> &dcorrection, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
