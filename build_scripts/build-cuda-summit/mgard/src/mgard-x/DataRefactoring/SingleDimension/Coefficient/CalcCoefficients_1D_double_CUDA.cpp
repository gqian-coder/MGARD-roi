/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/DataRefactoring/SingleDimension/Coefficient/CalcCoefficients.hpp"
// clang-format off
namespace mgard_x {

template void CalcCoefficients<1, double, CUDA>(
    DIM current_dim, SubArray<1, double, CUDA> ratio,
    SubArray<1, double, CUDA> v,
    SubArray<1, double, CUDA> coarse,
    SubArray<1, double, CUDA> coeff, int queue_idx);

} // namespace mgard_x
// clang-format on
