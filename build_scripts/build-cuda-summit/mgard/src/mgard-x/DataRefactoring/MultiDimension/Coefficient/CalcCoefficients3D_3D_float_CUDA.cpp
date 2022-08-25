/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Coefficient/CalcCoefficients3D.hpp"
// clang-format off
namespace mgard_x {

template void CalcCoefficients3D<3, float, CUDA>(
    Hierarchy<3, float, CUDA> &hierarchy,
    SubArray<3, float, CUDA> dinput,
    SubArray<3, float, CUDA> &doutput, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
