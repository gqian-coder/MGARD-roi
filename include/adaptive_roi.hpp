#ifndef ADAPTIVE_ROI
#define ADAPTIVE_ROI

#include <vector>
#include <cstddef>
#include <cstdint>

namespace mgard {

template <typename T> struct box_coord2d {
    T x0;
    T x1;
    T y0;
    T y1;
};

template <typename T> struct customized_hierarchy {
    T *level;
    T L;
    T Row;
    T Col;
    T Height;
    T l_th;
};

template <typename T>
void hist_blc_coord(std::vector<box_coord2d<T>> &blc_set, int nbin_R, int nbin_C,
                        int Row, int Col, size_t bin_w, size_t r0, size_t c0);


template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v);

template <typename T1, typename T2>
struct box_coord2d <T2> expand_boundary(struct box_coord2d<T2> coord, const T1* u_mc, size_t nb_w, T1 nb_thresh,
                                   customized_hierarchy <T2> &c_hierarchy);

template <typename T1, typename T2>
std::vector<struct box_coord2d<T2>> filter_hist_blc(const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy, std::vector<box_coord2d <T2>> &blc_set, const T1 thresh, size_t bin_w);


template <typename T1, typename T2>
void check_nearby_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                      int r, int c, int rad, std::vector<size_t> R2, u_char flag);

template <typename T1, typename T2>
void check_horizontal_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                        int r, int c, int rad, std::vector<size_t> R2, u_char flag);

template <typename T1, typename T2>
void check_vertical_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                        int r, int c, int rad, std::vector<size_t> R2, u_char flag);

template <typename T1, typename T2>
void check_edge_coeff(customized_hierarchy <T2> &c_hierarchy, T1 *u_map, int rmin,
                      int rmax, int cmin, int cmax, int max_rad, std::vector<size_t> R2);

template <typename T1, typename T2>
void dfs_amr(std::vector<struct box_coord2d<T2>> blc_set, const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy,
             const std::vector<T1> thresh, int &depth, size_t &bin_w, const size_t bin_min,
             const std::vector<size_t> ratio_bin, T1* u_map, const size_t max_rad, const std::vector<size_t> R2);

} // namespace mgard
#include "adaptive_roi.tpp"
#endif

