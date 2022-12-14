/*
 * Copyright 2022, Oak Ridge National Laboratory.
 * MGARD-X: MultiGrid Adaptive Reduction of Data Portable across GPUs and CPUs
 * Author: Jieyang Chen (chenj3@ornl.gov)
 * Date: March 17, 2022
 */

#include "mgard-x/CompressionHighLevel/CompressionHighLevel.hpp"
// clang-format off
namespace mgard_x {

template void decompress<CUDA>(const void *compressed_data,
                                         size_t compressed_size,
                                         void *&decompressed_data,
                                         Config config,
                                         bool output_pre_allocated);

template void decompress<CUDA>(const void *compressed_data,
                                         size_t compressed_size,
                                         void *&decompressed_data,
                                         bool output_pre_allocated);

template void
decompress<CUDA>(const void *compressed_data, size_t compressed_size,
                           void *&decompressed_data, data_type &dtype,
                           std::vector<SIZE> &shape, Config config,
                           bool output_pre_allocated);

template void
decompress<CUDA>(const void *compressed_data, size_t compressed_size,
                           void *&decompressed_data, data_type &dtype,
                           std::vector<SIZE> &shape, bool output_pre_allocated);

} // namespace mgard_x
// clang-format on
