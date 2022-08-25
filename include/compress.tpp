// Copyright 2017, Brown University, Providence, RI.
// MGARD: MultiGrid Adaptive Reduction of Data
// Authors: Mark Ainsworth, Ozan Tugluk, Ben Whitney
// Corresponding Author: Ben Whitney, Qing Liu
//
// See LICENSE for details.
#ifndef COMPRESS_TPP
#define COMPRESS_TPP

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>
#include <limits>
#include <numeric>
#include <type_traits>
#include <vector>

#include "MGARDConfig.hpp"
#include "TensorMultilevelCoefficientQuantizer.hpp"
#include "TensorNorms.hpp"
#include "compressors.hpp"
#include "decompose.hpp"
#include "format.hpp"
#include "quantize.hpp"
#include "shuffle.hpp"
#include "adaptive_roi.hpp"

namespace mgard {

using DEFAULT_INT_T = std::int64_t;

template <std::size_t N, typename Real>
CompressedDataset<N, Real>
compress(const TensorMeshHierarchy<N, Real> &hierarchy, Real *const v,
         const Real s, const Real tolerance) {
  const std::size_t ndof = hierarchy.ndof();
  Real *const u = new Real[ndof];
  shuffle(hierarchy, v, u);
  pb::Header header;
  populate_defaults(header);
  hierarchy.populate(header);
  decompose(hierarchy, header, u);
  {
    pb::ErrorControl &e = *header.mutable_error_control();
    e.set_mode(pb::ErrorControl::ABSOLUTE);
    if (s == std::numeric_limits<Real>::infinity()) {
      e.set_norm(pb::ErrorControl::L_INFINITY);
    } else {
      e.set_norm(pb::ErrorControl::S_NORM);
      e.set_s(s);
    }
    e.set_tolerance(tolerance);
  }
  MemoryBuffer<unsigned char> quantized = quantization_buffer(header, ndof);
  quantize(hierarchy, header, s, tolerance, u, quantized.data.get());
  MemoryBuffer<unsigned char> buffer =
      compress(header, quantized.data.get(), quantized.size);
  delete[] u;
  return CompressedDataset<N, Real>(hierarchy, header, s, tolerance,
                                    buffer.data.release(), buffer.size);
}

template <std::size_t N, typename Real>
DecompressedDataset<N, Real>
decompress(const CompressedDataset<N, Real> &compressed) {
  const std::size_t ndof = compressed.hierarchy.ndof();
  Real *const dequantized = new Real[ndof];
  Real *const v = new Real[ndof];
  MemoryBuffer<unsigned char> quantized =
      quantization_buffer(compressed.header, ndof);

  decompress(compressed.header, const_cast<void *>(compressed.data()),
             compressed.size(), quantized.data.get(), quantized.size);
  dequantize(compressed, quantized.data.get(), dequantized);

  recompose(compressed.hierarchy, compressed.header, dequantized);
  unshuffle(compressed.hierarchy, dequantized, v);
  delete[] dequantized;
  return DecompressedDataset<N, Real>(compressed, v);
}


template <std::size_t N, typename Real>
CompressedDataset<N, Real>
compress_roi(const TensorMeshHierarchy<N, Real> &hierarchy, Real *const v, const Real s,
         const Real tolerance, const std::vector<Real> thresh,
         const size_t bin_max, const std::vector<size_t> ratio_bin, const size_t l_th,
	 const char* filename, bool wr /*1 for write 0 for read*/) {
  const std::size_t ndof = hierarchy.ndof();
  Real *const u = new Real[ndof];
  shuffle(hierarchy, v, u);
  pb::Header header;
  populate_defaults(header);
  hierarchy.populate(header);
  decompose(hierarchy, header, u);

  // QG: create a map for adaptive compression
  Real *const unshuffled_u = static_cast<Real *>(std::malloc(ndof * sizeof(Real)));
  unshuffle(hierarchy, u, unshuffled_u);
  const std::array<std::size_t, N> &SHAPE = hierarchy.shapes.back();
  Real* u_map = static_cast<Real *>(std::malloc(ndof * sizeof(Real)));

  {
    pb::ErrorControl &e = *header.mutable_error_control();
    e.set_mode(pb::ErrorControl::ABSOLUTE);
    if (s == std::numeric_limits<Real>::infinity()) {
      e.set_norm(pb::ErrorControl::L_INFINITY);
    } else {
      e.set_norm(pb::ErrorControl::S_NORM);
      e.set_s(s);
    }
    e.set_tolerance(tolerance);
  }

  if ((wr==0) && (filename!=NULL)) {
    FILE *file = fopen(filename, "rb");
    fread(u_map, sizeof(Real), ndof, file);
  } else {
    int depth = 1;
    size_t bin_w = bin_max;
    size_t bin_min = bin_max;
    for (int i=0; i<ratio_bin.size(); i++) {
        bin_min = (size_t)std::ceil(bin_min / ratio_bin.at(i));
    }
    // std::memset(u_map, 0, ndof*sizeof(Real));
    //QG: Assume the array is 2D
    int Dim2=1, r, c, h;
    struct customized_hierarchy <size_t> c_hierarchy;
    c_hierarchy.level = new size_t[ndof];
    c_hierarchy.L     = hierarchy.L;
    c_hierarchy.Row   = (int)SHAPE[0];
    if (N > 3) {
        std::cout << "Adaptive compression does not support dim > 3!!!\n";
        return NULL; 
    }
    if (N==1) {
        c_hierarchy.Col    = 1;
        c_hierarchy.Height = 1;
    } else{
        c_hierarchy.Col = (int)SHAPE[1];
	if (N==3) {
            c_hierarchy.Height = (int)SHAPE[2];
            Dim2 = c_hierarchy.Height * c_hierarchy.Col;
	}
        else {
            c_hierarchy.Height = 1;
        }
    }
    // QG: get the level of each node (can be improved)
    for (std::size_t i=0; i<ndof; i++) {
        std::array<std::size_t, N> multiindex;
        if (N==1) {
            multiindex[0] = i;
        } else if (N==2) {
            r = int(i/c_hierarchy.Col);
            c = i - r*c_hierarchy.Col;
            multiindex[0] = r;
            multiindex[1] = c;
        } else {
            r = (int)(i / Dim2);
            c = (int)((i - r*Dim2) / c_hierarchy.Height);
            h = i % c_hierarchy.Height;
            multiindex[0] = r;
            multiindex[1] = c;
            multiindex[2] = h;
        }
        c_hierarchy.level[i] = hierarchy.date_of_birth(multiindex);
        // preserve coefficients w/ level < l_thresh w/ higher accuracy  
        if (c_hierarchy.level[i]<l_th) {
            u_map[i] = 125;
        } else {
            u_map[i] = 255;
        }
    }
    c_hierarchy.l_th = l_th;
    if (l_th<hierarchy.L) {
        // number of bins in the 1st layer of histogram
        int nbin_R = (int)std::ceil((float)c_hierarchy.Row / bin_max);
        int nbin_C = (N>=2) ? (int)std::ceil((float)c_hierarchy.Col / bin_max) : 1;
        int nbin_H = (N>=3) ? (int)std::ceil((float)c_hierarchy.Height / bin_max) : 1;
        std::vector<struct box_coord2d <size_t>> blc_set(nbin_R*nbin_C*nbin_H);
//        std::cout << "nbin_R = " << nbin_R << ", nbin_C = " << nbin_C << ", nbin_H = " << nbin_H << "\n";
//        std::cout << "Col = " << c_hierarchy.Col << ", Row = " << c_hierarchy.Row << ", Height = " << c_hierarchy.Height << "\n";
        size_t max_rad = (size_t)(2 * (1<<(c_hierarchy.L - c_hierarchy.l_th+1))); // 2nd peak to the 0th center
        std::vector<size_t>R2(c_hierarchy.L-c_hierarchy.l_th+1);
        R2.at(0) = max_rad;
//        std::cout << "Lmax = " << c_hierarchy.L << ",max_rad = " << max_rad << "\n";
        for (int i=1; i<c_hierarchy.L-c_hierarchy.l_th+1; i++) {
            R2.at(i) = R2.at(i-1) / 2;
        }
        hist_blc_coord<size_t>(blc_set, nbin_R, nbin_C, c_hierarchy.Row, c_hierarchy.Col, bin_max, 0, 0);
        std::vector<struct box_coord2d<size_t>>child_set = filter_hist_blc<Real>(unshuffled_u, c_hierarchy, blc_set, thresh[0], bin_max);

        // depth first search for hierachical block refinement
	dfs_amr<Real, size_t>(child_set, unshuffled_u, c_hierarchy, thresh, depth, bin_w, bin_min, ratio_bin, u_map, max_rad, R2);
//        std::cout << "number of blocks after the 1st thresholding: " << child_set.size() << "/" << blc_set.size() <<"\n";
    }
    if (filename != NULL) {
        FILE *fp = fopen (filename, "wb");
        fwrite (u_map , sizeof(Real), ndof, fp);
        fclose(fp);
    }
  }
  shuffle(hierarchy, u_map, unshuffled_u);

// QG: scalar of eb used for coefficients in non-RoI : RoI
// 2D case is bounded by the horizontal direction of coefficient_nodal error propagation
// 3rd peak, scalar=125 for 2D if nodal nodal and scalar=50 if nodal coefficient vertical
//  1D case: scalar = 22 for d=2h  (0.634 * w^2), w=0.268
//  2D case: scalar = 23 for d=2h  (0.5915 * w^2)
  size_t scalar = (N==1) ? 22 : ((N==2) ? 23 : 1);

  MemoryBuffer<unsigned char> quantized = quantization_buffer(header, ndof);
  quantize(hierarchy, header, s, tolerance, u, quantized.data.get());
  MemoryBuffer<unsigned char> buffer =
      compress(header, quantized.data.get(), quantized.size);
  delete[] u;
  return CompressedDataset<N, Real>(hierarchy, header, s, tolerance,
                                    buffer.data.release(), buffer.size);
}

} // namespace mgard
#endif
