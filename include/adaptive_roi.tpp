#include <algorithm>
#include <cmath>
#include "adaptive_roi.hpp"

namespace mgard {

template <typename T>
void hist_blc_coord(std::vector<box_coord2d <T>> &blc_set, int nbin_R, int nbin_C,
                        int Row, int Col, size_t bin_w, size_t r0, size_t c0) {
        int k=0;
        std::vector<size_t> coord_r0(nbin_R);
        std::vector<size_t> coord_r1(nbin_R);
        std::vector<size_t> coord_c0(nbin_C);
        std::vector<size_t> coord_c1(nbin_C);
        for (int r=0; r<nbin_R-1; r++) {
            coord_r0.at(r) = r0+r*bin_w;
            coord_r1.at(r) = coord_r0.at(r) + bin_w;
        }
        coord_r0.at(nbin_R-1) = r0+(nbin_R-1)*bin_w;
        coord_r1.at(nbin_R-1) = r0+Row;
    if (Col>1) { // 2D
        for (int c=0; c<nbin_C-1; c++) {
            coord_c0.at(c) = c0+c*bin_w;
            coord_c1.at(c) = coord_c0.at(c) + bin_w;
        }
        coord_c0.at(nbin_C-1) = c0+(nbin_C-1)*bin_w;
        coord_c1.at(nbin_C-1) = c0+Col;
    } else { // 1D
        std::fill(coord_c0.begin(), coord_c0.end(), 0);
        std::fill(coord_c1.begin(), coord_c1.end(), 1);
    }

    for(int r=0; r<nbin_R; r++) {
        for (int c=0; c<nbin_C; c++) {
            blc_set.at(k).x0 = coord_r0.at(r);
            blc_set.at(k).x1 = coord_r1.at(r);
            blc_set.at(k).y0 = coord_c0.at(c);
            blc_set.at(k).y1 = coord_c1.at(c);
            k++;
        }
    }
}

template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

//#define KAI 0.08
#define KAI 0
#define Alpha 1.414


template <typename T1, typename T2>
struct box_coord2d <T2> expand_boundary(struct box_coord2d <T2> coord, const T1* u_mc, size_t nb_w, T1 nb_thresh,
                                   customized_hierarchy <T2> &c_hierarchy) {
    size_t Row = c_hierarchy.Row;
    size_t Col = c_hierarchy.Col;
    struct box_coord2d <T2> child_blc = {coord.x0, coord.x1, coord.y0, coord.y1};
    int temp = coord.x0 - nb_w;
    size_t ext_r0 = (temp>0) ? temp : 0;
    temp = coord.x1 + nb_w;
    size_t ext_r1 = (temp<Row) ? temp : Row;
    temp = coord.y0 - nb_w;
    size_t ext_c0 = (temp>0) ? temp : 0;
    temp = coord.y1 + nb_w;
    size_t ext_c1 = (temp<Col) ? temp : Col;
    T1 nbw_top=0, nbw_bottom=0, nbw_left=0, nbw_right=0;
//    std::cout << "ext_r0=" << ext_r0 << ", ext_r1=" << ext_r1 << ", ext_c0=" << ext_c0 << ", ext_c1=" << ext_c1 << "\n";
//    std::cout << "{" << coord.x0 << ", " << coord.x1 << ", " << coord.y0 << ", " << coord.y1 << "}\n";
    size_t l;
    for (int r=ext_r0; r<ext_r1; r++) {
        int r0 = r*Col;
        if (r<coord.x0) {
            for (int c=r0+coord.y0; c<r0+coord.y1; c++){
                l = c_hierarchy.level[c];
                nbw_top += std::abs(u_mc[c]) / std::pow(Alpha, c_hierarchy.L-l);
            }
        } else if ((r>=coord.x0) && (r<coord.x1)) {
            for (int c=r0+ext_c0; c<r0+coord.y0; c++) {
                l = c_hierarchy.level[c];
                nbw_left  += std::abs(u_mc[c]) / std::pow(Alpha, c_hierarchy.L-l);
            }
            for (int c=r0+coord.y1; c<r0+ext_c1; c++) {
                l = c_hierarchy.level[c];
                nbw_right += std::abs(u_mc[c]) / std::pow(Alpha, c_hierarchy.L-l);
            }
        } else {
            for (int c=r0+coord.y0; c<r0+coord.y1; c++) {
                l = c_hierarchy.level[c];
                nbw_bottom += std::abs(u_mc[c]) / std::pow(Alpha, c_hierarchy.L-l);
            }
        }
    }
    if (nbw_top >= nb_thresh) {
        child_blc.x0 = ext_r0;
    }
    if (nbw_bottom >= nb_thresh) {
        child_blc.x1 = ext_r1;
    }
    if (nbw_left >= nb_thresh) {
        child_blc.y0 = ext_c0;
    }
    if (nbw_right >= nb_thresh) {
        child_blc.y1 = ext_c1;
    }
//      std::cout << "child_blc: {" << child_blc.x0 << ", " << child_blc.x1 << ", " << child_blc.y0 << ", " << child_blc.y1 << "}\n";
    return child_blc;
}

// return the id of blocks whose coefficients sum are large and expand the block boundary
template <typename T1, typename T2>
std::vector<struct box_coord2d <T2>> filter_hist_blc(const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy, std::vector<box_coord2d <T2>> &blc_set, const T1 thresh, size_t bin_w) {
    size_t nbins = blc_set.size();
//  std::cout << "number of bins: " << nbins << "\n";
    std::vector<T1>hist_w (nbins, 0);
    for (int i=0; i<nbins; i++) {
        for (int r=blc_set.at(i).x0; r<blc_set.at(i).x1; r++) {
            int r0 = r*c_hierarchy.Col;
            for (int c=blc_set.at(i).y0; c<blc_set.at(i).y1; c++) {
                size_t l = c_hierarchy.level[r0+c];
                hist_w[i] += std::abs(u_mc[r0+c]) / std::pow(Alpha, c_hierarchy.L-l);
            }
        }
        T1 area = (blc_set.at(i).x1 - blc_set.at(i).x0) * (blc_set.at(i).y1 - blc_set.at(i).y0);
        hist_w[i] = hist_w[i] / area;
    }
    std::vector<size_t> sid_umc = sort_indexes<T1>(hist_w);
//  for (int i=0; i<nbins; i++) {
//          std::cout << "u_mc[" << sid_umc[i] << "] = " << hist_w[sid_umc[i]] << "\n";
//  }
    int n_filtered = (int)std::ceil(thresh*nbins);
    std::vector<struct box_coord2d<T2>> filtered_set(n_filtered);
    // expand the boundary by KAI
    size_t nb_w = (size_t)(std::ceil((float)bin_w * KAI));
    float edge_ratio = (float)nb_w / (float)bin_w;
//  std::cout << "bin_w = " << bin_w << ", edge ratio: " << edge_ratio << "\n";
    for (int i=0; i<n_filtered; i++) {
        if ((bin_w > 1) && (KAI>0)) {
            T1 nb_thresh = edge_ratio*hist_w.at(sid_umc.at(i));
            filtered_set.at(i) = expand_boundary<T1, T2>(blc_set.at(sid_umc.at(i)), u_mc, nb_w, nb_thresh, c_hierarchy);
        } else {
            filtered_set.at(i) = blc_set.at(sid_umc.at(i));
        }
    }
    return filtered_set;
}

template <typename T1, typename T2>
void check_nearby_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                      int r, int c, int rad, std::vector<size_t> R2, u_char flag)
{
  int k, temp;
  std::vector<int> top_r{0, 0, r, r};
  std::vector<int> bottom_r{r+1, r+1, c_hierarchy.Row, c_hierarchy.Row};
  std::vector<int> left_c{0, c, 0, c};
  std::vector<int> right_c{c+1, c_hierarchy.Col, c+1, c_hierarchy.Col};
  size_t l, backoff_dist;
  temp = r-rad;
  if (temp>0) {top_r[0] = temp; top_r[1] = temp;}
  temp = c-rad;
  if (temp>0) {left_c[0] = temp; left_c[2] = temp;}
  temp = r+rad+1;
  if (temp<c_hierarchy.Row) {bottom_r[2] = temp; bottom_r[3] = temp;}
  temp = c+rad+1;
  if (temp<c_hierarchy.Col) {right_c[1] = temp; right_c[3] = temp;}
  size_t l_th = c_hierarchy.l_th;
  size_t Col  = c_hierarchy.Col;
  int check_n = 0;
  for (int fg=0; fg<4; fg++) { // top left, top right, bottom left, bottom right
    if ((flag & (1<<fg))) {
      for (int rr=top_r[fg]; rr<bottom_r[fg]; rr++) {
        for (int cc=left_c[fg]; cc<right_c[fg]; cc++) {
          k = rr*Col+cc;
          if (u_map[k] !=0) {
            l = c_hierarchy.level[k];
            if (l>=l_th) {
              backoff_dist = R2[l-l_th];
              if (std::abs(rr-r)+std::abs(cc-c) < backoff_dist) {
                u_map[k] = 125;
              }
            }
          }
        }
      }
    }
  }
}

template <typename T1, typename T2>
void check_horizontal_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                        int r, int c, int rad, std::vector<size_t> R2, u_char flag)
{
    int k;
    size_t l, backoff_dist;
    size_t l_th = c_hierarchy.l_th;
    size_t Col  = c_hierarchy.Col;
    int cmin = flag ? c : ((c-rad>0) ? (c-rad):0);;//(c+5) : ((c-rad>0) ? (c-rad):0);
    int cmax = flag ? ((c+rad<Col) ? (c+rad):Col-1) : c;// (c-5);
    for (int cc=cmin+1; cc<=cmax; cc++) {
        k = r*Col + cc;
        if (u_map[k] == 255) { // l<l_th already being set to buffer zone
            l = c_hierarchy.level[k];
            backoff_dist = R2[l-l_th];
            if (std::abs(cc-c) < backoff_dist) {
                u_map[k] = 125;
            }
        }
    }
}

template <typename T1, typename T2>
void check_vertical_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                        int r, int c, int rad, std::vector<size_t> R2, u_char flag)
{
    int k;
    int rmin = flag ? r : ((r-rad>0) ? (r-rad):0); //(r+5) : ((r-rad>0) ? (r-rad):0);
    int rmax = flag ? ((r+rad<c_hierarchy.Row) ? (r+rad):c_hierarchy.Row-1) : c;//(r-5);
    size_t l, backoff_dist;
    size_t l_th = c_hierarchy.l_th;
    size_t Col  = c_hierarchy.Col;
    for (int rr=rmin+1; rr<=rmax; rr++) {
        k = rr*Col + c;
        if (u_map[k] == 255) { // l<l_th already being set to buffer zone
            l = c_hierarchy.level[k];
            backoff_dist = R2[l-l_th];
            if (std::abs(rr-r) < backoff_dist) {
                u_map[k] = 125;
            }
        }
    }
}

//  expand the 2D region surrounding the selected coefficients
template <typename T1, typename T2>
void check_edge_coeff(customized_hierarchy <T2> &c_hierarchy, T1 *u_map, int rmin,
                      int rmax, int cmin, int cmax, int max_rad, std::vector<size_t> R2)
{
    std::vector<u_char>flag {1,2,4,8};//{7, 11, 13, 14};
    // the top and bottom horizontal edge
    size_t k, Col, Row, rstart;
    Col = c_hierarchy.Col;
    Row = c_hierarchy.Row;
    rstart = ((size_t)rmin) * Col;
    for (int cc=cmin+1; cc<cmax; cc++) {
        k = rstart + cc;
        if (u_map[k]!=0)
            check_vertical_lgc<T1>(c_hierarchy, u_map, rmin, cc, max_rad, R2, 0);
    }
    if (u_map[rstart+cmin]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, rmin, cmin, max_rad, R2, flag[0]);

    if (u_map[rstart+cmax]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, rmin, cmax, max_rad, R2, flag[1]);
    rstart = ((size_t)rmax) * Col;
    for (int cc=cmin+1; cc<cmax; cc++) {
        k = rstart + cc;
        if (u_map[k]!=0)
            check_vertical_lgc<T1, T2>(c_hierarchy, u_map, rmax, cc, max_rad, R2, 1);
    }
    if (u_map[rstart+cmin]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, rmax, cmin, max_rad, R2, flag[2]);
    if (u_map[rstart+cmax]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, rmax, cmax, max_rad, R2, flag[3]);

    // the left and right vertical edge
    for (int rr=rmin+1; rr<rmax; rr++) {
        k = ((size_t)rr*Col + cmin);
        if (u_map[k]!=0)
            check_horizontal_lgc<T1, T2>(c_hierarchy, u_map, rr, cmin, max_rad, R2, 0);
    }
    for (int rr=rmin+1; rr<rmax; rr++) {
        k = ((size_t)rr*Col + cmax);
        if (u_map[k]!=0)
            check_horizontal_lgc<T1, T2>(c_hierarchy, u_map, rr, cmax, max_rad, R2, 1);
    }
}


template <typename T1, typename T2>
void dfs_amr(std::vector<struct box_coord2d <T2>> blc_set, const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy,
             const std::vector<T1> thresh, int &depth, size_t &bin_w, const size_t bin_min,
             const std::vector<size_t> ratio_bin, T2* u_map, const size_t max_rad, const std::vector<size_t> R2) {
//      std::cout << "initial bin_w = " << bin_w << ", ratio_R = " << ratio_bin[depth] << "\n";
    bin_w = std::ceil(bin_w / ratio_bin[depth]);
//  std::cout << "depth = " << depth << ", bin_w = " << bin_w << "\n";
//  std::cout << "number of bins: " << blc_set.size() << "\n";
    for (int k=0; k<blc_set.size(); k++) {
        size_t bR  = blc_set.at(k).x1 - blc_set.at(k).x0;
        size_t bC  = blc_set.at(k).y1 - blc_set.at(k).y0;
        size_t nbin_R = (size_t)std::ceil((float)bR / bin_w);
        size_t nbin_C = (size_t)std::ceil((float)bC / bin_w);
//      std::cout << "{" <<  blc_set.at(k).x0 << ", " <<  blc_set.at(k).x1 << ", " <<  blc_set.at(k).y0 << ", " <<  blc_set.at(k).y1 << "}\n";
//      std::cout << "k = " << k << ", bR = " << bR << ", bC = " << bC << ", nbin_R = " << nbin_R << ", nbin_C = " << nbin_C << "\n";
        std::vector<struct box_coord2d<T2>> child_set(nbin_R*nbin_C);
        hist_blc_coord<T2>(child_set, nbin_R, nbin_C, bR, bC, bin_w, blc_set.at(k).x0, blc_set.at(k).y0);
        std::vector<struct box_coord2d<T2>>filtered_set = filter_hist_blc<T1, T2>(u_mc, c_hierarchy, child_set, thresh[depth], bin_w);
        if (bin_w <= bin_min){
//          std::cout << "max_rad: " << max_rad << "\n";
            for (int cid=0; cid<filtered_set.size(); cid++) {
                int r, c;
                for (r=filtered_set.at(cid).x0; r<filtered_set.at(cid).x1; r++) {
                    int r0 = r * c_hierarchy.Col;
                    for (c=filtered_set.at(cid).y0; c<filtered_set.at(cid).y1; c++) {
                        std::size_t l = c_hierarchy.level[r0+c];
                        if (l >= c_hierarchy.l_th) {
                            int rad = static_cast<int>(1<<(c_hierarchy.L - l));
                            int rmin  = (r-rad > 0) ? (r-rad) : 0;
                            int rmax = (r+rad+1<c_hierarchy.Row) ? (r+rad+1) : c_hierarchy.Row;
                            int cmin  = (c-rad>0) ? c-rad : 0;
                            int cmax = (c+rad+1<c_hierarchy.Col) ? (c+rad+1) : c_hierarchy.Col;
                            check_edge_coeff(c_hierarchy, u_map, rmin, rmax, cmin, cmax, max_rad, R2);
                            for (int rr=rmin; rr<rmax; rr++) {
                                int r00 = rr*c_hierarchy.Col;
                                for (int cc=cmin; cc<cmax; cc++) {
                                    u_map[r00+cc] = 0;
                                }
                            }
                        }
                    }
                }
            }
        } else {
//         std::cout << "call recursive dfs_arm, depth = " << depth << ", bin_w = " << bin_w << ", # of blocks: " << filtered_set.size() << " / " << child_set.size() << "\n";
            depth ++;
            dfs_amr<T1, T2>(filtered_set, u_mc, c_hierarchy, thresh, depth, bin_w, bin_min, ratio_bin, u_map, max_rad, R2);
        }
    }
    bin_w = size_t(bin_w * ratio_bin[depth]);
    depth --;
}

} // namespace mgard 
