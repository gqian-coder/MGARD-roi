/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/Coefficient/CalcCoefficientsND.hpp"
// clang-format off
namespace mgard_x {

template void CalcCoefficientsND<3, double, SERIAL>(
    Hierarchy<3, double, SERIAL> &hierarchy,
    SubArray<3, double, SERIAL> dinput1,
    SubArray<3, double, SERIAL> dinput2,
    SubArray<3, double, SERIAL> &doutput, SIZE l,
    int queue_idx);

} // namespace mgard_x
// clang-format on
