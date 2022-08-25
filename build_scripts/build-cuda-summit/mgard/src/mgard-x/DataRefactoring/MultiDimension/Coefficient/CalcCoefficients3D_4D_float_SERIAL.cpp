/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Coefficient/CalcCoefficients3D.hpp"
// clang-format off
namespace mgard_x {

template void CalcCoefficients3D<4, float, SERIAL>(
    Hierarchy<4, float, SERIAL> &hierarchy,
    SubArray<4, float, SERIAL> dinput,
    SubArray<4, float, SERIAL> &doutput, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
