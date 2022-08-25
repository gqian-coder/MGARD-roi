/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/CompressionLowLevel/CompressionLowLevel.hpp"
// clang-format off
namespace mgard_x {

template Array<1, unsigned char, CUDA>
compress<3, float, CUDA>(
    Hierarchy<3, float, CUDA> &hierarchy,
    Array<3, float, CUDA> &in_array,
    enum error_bound_type type, float tol, float s,
    float &norm, Config config);

} // namespace mgard_x
// clang-format on
