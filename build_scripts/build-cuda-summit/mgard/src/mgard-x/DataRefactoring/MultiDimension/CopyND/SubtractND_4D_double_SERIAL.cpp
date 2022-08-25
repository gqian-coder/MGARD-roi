/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/MultiDimension/CopyND/SubtractND.hpp"
// clang-format off
namespace mgard_x {

template void SubtractND<4, double, SERIAL>(
    SubArray<4, double, SERIAL> dinput,
    SubArray<4, double, SERIAL> &doutput, int queue_idx);

} // namespace mgard_x
// clang-format on
