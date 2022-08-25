/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Coefficient/CalcCoefficients3D.hpp"
// clang-format off
namespace mgard_x {

template void CalcCoefficients3D<1, double, CUDA>(
    Hierarchy<1, double, CUDA> &hierarchy,
    SubArray<1, double, CUDA> dinput,
    SubArray<1, double, CUDA> &doutput, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
