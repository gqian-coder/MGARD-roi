/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/Coefficient/CoefficientsRestore.hpp"
// clang-format off
namespace mgard_x {

template void CoefficientsRestore<1, double, SERIAL>(
    DIM current_dim, SubArray<1, double, SERIAL> ratio,
    SubArray<1, double, SERIAL> v,
    SubArray<1, double, SERIAL> coarse,
    SubArray<1, double, SERIAL> coeff, int queue_idx);

} // namespace mgard_x
// clang-format on
