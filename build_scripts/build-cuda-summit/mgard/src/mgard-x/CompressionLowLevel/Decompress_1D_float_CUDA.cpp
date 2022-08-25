/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/CompressionLowLevel/CompressionLowLevel.hpp"
// clang-format off
namespace mgard_x {

template Array<1, float, CUDA>
decompress<1, float, CUDA>(
    Hierarchy<1, float, CUDA> &hierarchy,
    Array<1, unsigned char, CUDA> &compressed_array,
    enum error_bound_type type, float tol, float s,
    float norm, Config config);

} // namespace mgard_x
// clang-format on
