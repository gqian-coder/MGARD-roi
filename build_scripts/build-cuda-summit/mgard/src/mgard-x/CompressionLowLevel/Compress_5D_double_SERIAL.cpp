/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/CompressionLowLevel/CompressionLowLevel.hpp"
// clang-format off
namespace mgard_x {

template Array<1, unsigned char, SERIAL>
compress<5, double, SERIAL>(
    Hierarchy<5, double, SERIAL> &hierarchy,
    Array<5, double, SERIAL> &in_array,
    enum error_bound_type type, double tol, double s,
    double &norm, Config config);

} // namespace mgard_x
// clang-format on
