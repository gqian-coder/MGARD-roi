/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/CompressionHighLevel/CompressionHighLevel.hpp"
// clang-format off
namespace mgard_x {

template void compress<SERIAL>(DIM D, data_type dtype,
                                       std::vector<SIZE> shape, double tol,
                                       double s, enum error_bound_type mode,
                                       const void *original_data,
                                       void *&compressed_data,
                                       size_t &compressed_size, Config config,
                                       bool output_pre_allocated);

template void
compress<SERIAL>(DIM D, data_type dtype, std::vector<SIZE> shape,
                         double tol, double s, enum error_bound_type mode,
                         const void *original_data, void *&compressed_data,
                         size_t &compressed_size, bool output_pre_allocated);

template void compress<SERIAL>(
    DIM D, data_type dtype, std::vector<SIZE> shape, double tol, double s,
    enum error_bound_type mode, const void *original_data,
    void *&compressed_data, size_t &compressed_size,
    std::vector<const Byte *> coords, Config config, bool output_pre_allocated);

template void compress<SERIAL>(
    DIM D, data_type dtype, std::vector<SIZE> shape, double tol, double s,
    enum error_bound_type mode, const void *original_data,
    void *&compressed_data, size_t &compressed_size,
    std::vector<const Byte *> coords, bool output_pre_allocated);
} // namespace mgard_x
// clang-format on
