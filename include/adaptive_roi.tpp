#include <algorithm>
#include <cmath>
#include "adaptive_roi.hpp"

namespace mgard {

template <typename T>
void hist_blc_coord(std::vector<cube_<T>> &blc_set /*top left coord*/, 
                    cube_<T> nbin, cube_<T> bin_w, cube_<T>init) 
{
    size_t k, r, c, h, dim2;
    dim2 = nbin.c*nbin.h;
    for (r=0; r<nbin.r; r++) {
        for (c=0; c<nbin.c; c++) {
            for (h=0; h<nbin.h; h++) {
                k = r*dim2 + c*nbin.h + h;
                blc_set.at(k).r = init.r + r * bin_w.r; 
                blc_set.at(k).c = init.c + c * bin_w.c;
                blc_set.at(k).h = init.h + h * bin_w.h;
            }
        }
    }
}

template <typename T>
size_t blc_coord_gb(std::vector<cube_<int>> &c_blc, std::vector<cube_<int>>p_blc, 
                    size_t nblc, std::vector<cube_<T>> bin_w, size_t depth)
{
    size_t k=0, r, c, h;
    struct cube_<size_t> nbin;
    T row = bin_w[0].r;
    T col = bin_w[0].c;
    T hig = bin_w[0].h;
    struct cube_<T> prev_bw = {bin_w[depth-1].r, bin_w[depth-1].c, bin_w[depth-1].h};
    struct cube_<T> curr_bw = {bin_w[depth].r, bin_w[depth].c, bin_w[depth].h};
    for (size_t i=0; i<nblc; i++) {
        nbin.r = (p_blc[i].r+prev_bw.r>row) ? (size_t)std::ceil(((float)(row-p_blc[i].r))/curr_bw.r) : (size_t)std::ceil(((float)prev_bw.r)/curr_bw.r);
        nbin.c = (p_blc[i].c+prev_bw.c>col) ? (size_t)std::ceil(((float)(col-p_blc[i].c))/curr_bw.c) : (size_t)std::ceil(((float)prev_bw.c)/curr_bw.c);
        nbin.h = (p_blc[i].h+prev_bw.h>hig) ? (size_t)std::ceil(((float)(hig-p_blc[i].h))/curr_bw.h) : (size_t)std::ceil(((float)prev_bw.h)/curr_bw.h);
        for (r=0; r<nbin.r; r++) {
            for (c=0; c<nbin.c; c++) {
                for (h=0; h<nbin.h; h++) {
                    c_blc.at(k).r = p_blc[i].r + r * curr_bw.r;
                    c_blc.at(k).c = p_blc[i].c + c * curr_bw.c;
                    c_blc.at(k).h = p_blc[i].h + h * curr_bw.h;               
                    k++;
                }
            }
        }
    } 
    return k;
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

#define Alpha 1.414

// return the id of blocks whose coefficients sum are large and expand the block boundary
template <typename T1, typename T2>
void filter_hist_blc(const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy, std::vector<cube_<int>> &blc_set, 
                    std::vector<cube_<int>> &filtered_set, const T1 thresh, cube_<T2> bin_w, size_t &nbins)
{
//  std::cout << "number of bins: " << nbins << "\n";
    std::vector<T1>hist_w (nbins, 0);
    size_t dim2 = c_hierarchy.Col * c_hierarchy.Height;
    size_t k;
    T1 area, r1, c1, h1;
    for (size_t i=0; i<nbins; i++) {
        r1 = blc_set.at(i).r + bin_w.r;
        area = bin_w.r;
        if (r1>c_hierarchy.Row) { 
            area = r1-c_hierarchy.Row;
            r1   = c_hierarchy.Row;
        }
        c1 = blc_set.at(i).c + bin_w.c;
        if (c1>c_hierarchy.Col) {
            area = area * (c1 - c_hierarchy.Col);
            c1   = c_hierarchy.Col;
        } else {
            area = area * bin_w.c; 
        }
        h1 = blc_set.at(i).h + bin_w.h;
        if (h1>c_hierarchy.Height) { 
            area = area * (h1 - c_hierarchy.Height);
            h1   = c_hierarchy.Height;
        } else {
            area = area * bin_w.h;
        }
        for (T2 r=blc_set.at(i).r; r<r1; r++) {
            T2 r0 = r*dim2;
            for (T2 c=blc_set.at(i).c; c<c1; c++) {
                T2 c0 = c*c_hierarchy.Height;
                for (T2 h=blc_set.at(i).h; h<h1; h++) {
                    k = r0+c0+h;
                    T2 l = c_hierarchy.level[k];
                    if (l<c_hierarchy.l_th)
                        hist_w[i] += std::abs(u_mc[k]) / std::pow(Alpha, c_hierarchy.L-l);
                }
            }
        } 
        // normalize to remove the edge effect
        hist_w[i] = hist_w[i] / area;
    }
    std::vector<size_t> sid_umc = sort_indexes<T1>(hist_w);
    nbins = (size_t)std::ceil(thresh*nbins); 
    for (int i=0; i<nbins; i++) {
            filtered_set.at(i) = blc_set.at(sid_umc.at(i));
    }
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
                u_map[k] = BUFFER_ZONE;
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
        if (u_map[k] == BACKGROUND) { // l<l_th already being set to buffer zone
            l = c_hierarchy.level[k];
            backoff_dist = R2[l-l_th];
            if (std::abs(cc-c) < backoff_dist) {
                u_map[k] = BUFFER_ZONE;
            }
        }
    }
}

template <typename T1, typename T2>
void check_vertical_lgc(customized_hierarchy <T2> &c_hierarchy, T1 *u_map,
                        int r, int c, int rad, std::vector<size_t> R2, u_char flag)
{
    int k;
    int rmin = flag ? r : ((r-rad>0) ? (r-rad):0); 
    int rmax = flag ? ((r+rad<c_hierarchy.Row) ? (r+rad):c_hierarchy.Row-1) : r;
    size_t l, backoff_dist;
    size_t l_th = c_hierarchy.l_th;
    size_t Col  = c_hierarchy.Col;
    for (int rr=rmin+1; rr<=rmax; rr++) {
        k = rr*Col + c;
        if (u_map[k] == BACKGROUND) { // l<l_th already being set to buffer zone
            l = c_hierarchy.level[k];
            backoff_dist = R2[l-l_th];
            if (std::abs(rr-r) < backoff_dist) {
                u_map[k] = BUFFER_ZONE;
            }
        }
    }
}

//  expand the 2D region surrounding the selected coefficients
template <typename T1, typename T2>
void check_edge_coeff(customized_hierarchy <T2> &c_hierarchy, T1 *u_map, 
                      struct cube_<T2> blc_min, struct cube_<T2> blc_max,
                      std::vector<size_t> R2)
{
    std::vector<u_char>flag {1,2,4,8};
    // the top and bottom horizontal edge
    size_t k, Col, Row, rstart;
    Col = c_hierarchy.Col;
    Row = c_hierarchy.Row;
    rstart = ((size_t)blc_min.r) * Col;
    for (int cc=blc_min.c+1; cc<blc_max.c; cc++) {
        k = rstart + cc;
        if (u_map[k]!=0)
            check_vertical_lgc<T1, T2>(c_hierarchy, u_map, blc_min.r, cc, R2[0], R2, 0);
    }
    if (u_map[rstart+blc_min.c]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, blc_min.r, blc_min.c, R2[0], R2, flag[0]);

    if (u_map[rstart+blc_max.c]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, blc_min.r, blc_max.c, R2[0], R2, flag[1]);
    rstart = ((size_t)blc_max.r) * Col;
    for (int cc=blc_min.c+1; cc<blc_max.c; cc++) {
        k = rstart + cc;
        if (u_map[k]!=0)
            check_vertical_lgc<T1, T2>(c_hierarchy, u_map, blc_max.r, cc, R2[0], R2, 1);
    }
    if (u_map[rstart+blc_min.c]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, blc_max.r, blc_min.c, R2[0], R2, flag[2]);
    if (u_map[rstart+blc_max.c]!=0)
        check_nearby_lgc<T1, T2>(c_hierarchy, u_map, blc_max.r, blc_max.c, R2[0], R2, flag[3]);

    // the left and right vertical edge
    for (int rr=blc_min.r+1; rr<blc_max.r; rr++) {
        k = ((size_t)rr*Col + blc_min.c);
        if (u_map[k]!=0)
            check_horizontal_lgc<T1, T2>(c_hierarchy, u_map, rr, blc_min.c, R2[0], R2, 0);
    }
    for (int rr=blc_min.r+1; rr<blc_max.r; rr++) {
        k = ((size_t)rr*Col + blc_max.c);
        if (u_map[k]!=0)
            check_horizontal_lgc<T1, T2>(c_hierarchy, u_map, rr, blc_max.c, R2[0], R2, 1);
    }
}


template <typename T1, typename T2>
void dfs_amr_2d(std::vector<struct cube_ <T2>> blc_set, const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy, 
            const std::vector<T1> thresh, int &depth, std::vector<cube_<T2>> bin_w, T1* u_map, 
            const std::vector<size_t> R2) 
{
//  std::cout << "number of bins: " << blc_set.size() << "\n";
    for (int k=0; k<blc_set.size(); k++) {
        struct cube_<T2> nbin;
        struct cube_<T2> prev_bin = bin_w[depth-1];
        struct cube_<T2> init_pos = {blc_set.at(k).r, blc_set.at(k).c, blc_set.at(k).h};
        if (blc_set.at(k).r+prev_bin.r>c_hierarchy.Row) prev_bin.r = c_hierarchy.Row-blc_set.at(k).r;
        if (blc_set.at(k).c+prev_bin.r>c_hierarchy.Col) prev_bin.c = c_hierarchy.Col-blc_set.at(k).c;
        if (blc_set.at(k).h+prev_bin.r>c_hierarchy.Height) prev_bin.h = c_hierarchy.Height-blc_set.at(k).h;
        nbin.r = (T2) std::ceil((float)prev_bin.r / bin_w[depth].r);
        nbin.c = (T2) std::ceil((float)prev_bin.c / bin_w[depth].c);
        nbin.h = (T2) std::ceil((float)prev_bin.h / bin_w[depth].h);
        size_t nblc = nbin.r*nbin.c*nbin.h;
        std::vector<struct cube_<int>> child_set(nblc);
        hist_blc_coord<T2>(child_set, nbin, bin_w[depth], init_pos);
        int n_filtered = (int)std::ceil(thresh*nblc);
        std::vector<struct cube_<int>> filtered_set(n_filtered);
        filter_hist_blc<T1, T2>(u_mc, c_hierarchy, child_set, filtered_set, thresh[depth], bin_w[depth], nblc);
        if (depth==bin_w.size()-1){
            for (int cid=0; cid<nblc; cid++) {
                T2 r, c, h;
                T2 r1 = filtered_set.at(cid).r + bin_w[depth].r;
                if (r1>c_hierarchy.Row) r1 = c_hierarchy.Row;
                T2 c1 = filtered_set.at(cid).c + bin_w[depth].c;
                if (c1>c_hierarchy.Col) c1 = c_hierarchy.Col;
                T2 h1 = filtered_set.at(cid).h + bin_w[depth].h;
                if (h1>c_hierarchy.Height) h1 = c_hierarchy.Height;
                T2 dim2 = c_hierarchy.Col * c_hierarchy.Height;
                for (r=filtered_set.at(cid).r; r<r1; r++) {
                    T2 r0 = r * dim2;
                    for (c=filtered_set.at(cid).c; c<c1; c++) {
                        T2 c0 = c * c_hierarchy.Height;
                        for (h=filtered_set.at(cid).h; h<h1; h++) {
                            T2 k = r0 + c0 + h;
                            std::size_t l = c_hierarchy.level[k];
                            if ((l >= c_hierarchy.l_th) && (u_map[k]>0)) {
                                T2 rad  = static_cast<int>(1<<(c_hierarchy.L - l));
                                struct cube_<T2> blc_min, blc_max;
                                blc_min.r = (r-rad > 0) ? (r-rad) : 0;
                                blc_max.r = (r+rad+1<c_hierarchy.Row) ? (r+rad+1) : c_hierarchy.Row;
                                blc_min.c = (c-rad>0) ? c-rad : 0;
                                blc_max.c = (c+rad+1<c_hierarchy.Col) ? (c+rad+1) : c_hierarchy.Col;
                                blc_min.h = (h-rad>0) ? h-rad : 0;
                                blc_max.h = (h+rad+1<c_hierarchy.Height) ? (h+rad+1) : c_hierarchy.Height;
                                check_edge_coeff(c_hierarchy, u_map, blc_min, blc_max, R2);
                                for (T2 rr=blc_min.r; rr<blc_max.r; rr++) {
                                    T2 r00 = rr*dim2;
                                    for (T2 cc=blc_min.c; cc<blc_max.c; cc++) {
                                        T2 c00 = cc*c_hierarchy.Height;
                                        for (T2 hh=blc_min.h; hh<blc_max.h; hh++) {
                                            u_map[r00+c00+hh] = 0;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
//         std::cout << "call recursive dfs_arm, depth = " << depth << ", bin_w = " << bin_w << ", # of blocks: " << filtered_set.size() << " / " << child_set.size() << "\n";
            depth ++;
            dfs_amr_2d<T1, T2>(filtered_set, u_mc, c_hierarchy, thresh, depth, bin_w, u_map, R2);
        }
    }
    depth --;
}

template <typename T1, typename T2>
bool buffer_zone(customized_hierarchy <T2> &c_hierarchy, T1* u_map, struct cube_<int> pos,
                const int rad)
{
    int rr, cc, hh;
    size_t r00, c00;
    size_t delta_r, delta_c, delta_h;
    struct cube_<size_t> epa_min, epa_max;
    size_t dim2 = c_hierarchy.Col*c_hierarchy.Height;
    epa_min.r = (pos.r-rad>0) ? (pos.r-rad) : 0;
    epa_max.r = (pos.r+rad+1<c_hierarchy.Row) ? (pos.r+rad+1) : c_hierarchy.Row;
    epa_min.c = (pos.c-rad>0) ? (pos.c-rad) : 0;
    epa_max.c = (pos.c+rad+1<c_hierarchy.Col) ? (pos.c+rad+1) : c_hierarchy.Col;
    epa_min.h = (pos.h-rad>0) ? (pos.h-rad) : 0;
    epa_max.h = (pos.h+rad+1<c_hierarchy.Height) ? (pos.h+rad+1) : c_hierarchy.Height;
    for (rr=epa_min.r; rr<epa_max.r; rr++) {
        r00 = rr*dim2;
        delta_r = std::abs(rr - pos.r);
        for (cc=epa_min.c; cc<epa_max.c; cc++) {
            c00 = cc*c_hierarchy.Height;
            delta_c = std::abs(cc - pos.c);
            for (hh=epa_min.h; hh<epa_max.h; hh++) {
                delta_h = std::abs(hh - pos.h);
                if ((delta_r+delta_c+delta_h<=rad) && (u_map[r00+c00+hh]==ROI)) 
                    return true;
            }
        }
    }
    return false;
}

// bin_w: [{Row, Col, Height}, ...], size == depth+1
//        minimally bin_w will have two levels
// depth: thresh.size()
template <typename T1, typename T2>
void amr_gb(const T1 *u_mc, customized_hierarchy <T2> &c_hierarchy, const std::vector<T1> thresh, 
            const std::vector<cube_<T2>> bin_w, T1* u_map, const std::vector<size_t> R2)
{
    size_t depth = thresh.size(); 
    // initialize parent blocks (maximal #)
    size_t nr = (size_t)std::ceil((float)bin_w[0].r / bin_w[depth].r);
    size_t nc = (size_t)std::ceil((float)bin_w[0].c / bin_w[depth].c);
    size_t nh = (size_t)std::ceil((float)bin_w[0].h / bin_w[depth].h);
    // initialize children blocks (maximal #)
    nr = (size_t)std::ceil((float)bin_w[0].r / bin_w[depth].r);
    nc = (size_t)std::ceil((float)bin_w[0].c / bin_w[depth].c);
    nh = (size_t)std::ceil((float)bin_w[0].h / bin_w[depth].h);
    std::vector<struct cube_<int>> child_blc(nr*nc*nh);
    std::vector<struct cube_<int>> parent_blc(nr*nc*nh);
    
    // mesh refinement
    size_t n_pblc = 1;
    struct cube_<int> blc_min, blc_max, epa_min, epa_max;
    for (size_t d=0; d<depth; d++) {
        n_pblc = blc_coord_gb<T2>(child_blc, parent_blc, n_pblc, bin_w, d+1);
        filter_hist_blc<T1, T2>(u_mc, c_hierarchy, child_blc, parent_blc, thresh[d], bin_w[d+1], n_pblc);
    } // parent_blc contains the filtered blocks
//    std::cout << "final bins..." << n_pblc << "\n";
    
    // RoI expansion  
    size_t dim2 = (size_t)c_hierarchy.Col * c_hierarchy.Height;
    size_t r0, c0, r00, c00,  k, l, epa_h;
    int rad, rr, cc, hh;
    for (size_t i=0; i<n_pblc; i++) {
        blc_min = parent_blc[i];
        blc_max.r = blc_min.r + bin_w[depth].r;
        blc_max.c = blc_min.c + bin_w[depth].c;
        blc_max.h = blc_min.h + bin_w[depth].h; 
//        std::cout << "blc_min.r = " << blc_min.r << ", .c = " << blc_min.c << ", .h = " << blc_min.h << "\n";
//        std::cout << "blc_max.r = " << blc_max.r << ", .c = " << blc_max.c << ", .h = " << blc_max.h << "\n";
        for (int r=blc_min.r; r<blc_max.r; r++) {
            r0 = r*dim2; 
            for (int c=blc_min.c; c<blc_max.c; c++) {
                c0 = c*c_hierarchy.Height;
                for (int h=blc_min.h; h<blc_max.h; h++) {
                    k = r0 + c0 + h;
                    l = c_hierarchy.level[k];
                    if (l>=c_hierarchy.l_th) { 
                        rad  = static_cast<int>(1<<(c_hierarchy.L - l));
                        epa_min.r = (r-rad > 0) ? (r-rad) : 0;
                        epa_max.r = (r+rad+1<c_hierarchy.Row) ? (r+rad+1) : c_hierarchy.Row;
                        epa_min.c = (c-rad>0) ? c-rad : 0;
                        epa_max.c = (c+rad+1<c_hierarchy.Col) ? (c+rad+1) : c_hierarchy.Col;
                        epa_min.h = (h-rad>0) ? h-rad : 0;
                        epa_max.h = (h+rad+1<c_hierarchy.Height) ? h+rad+1 : c_hierarchy.Height;
                        epa_h     = epa_max.h - epa_min.h;
                        for (rr=epa_min.r; rr<epa_max.r; rr++) {
                            r00 = rr*dim2;
                            for (cc=epa_min.c; cc<epa_max.c; cc++) {
                                c00 = cc*c_hierarchy.Height;
                                memset(&u_map[r00+c00+epa_min.h], 0, epa_h * sizeof(T1));
                            }
                        }
                    } 
                }
            }
        } 
    }
//    std::cout << "RoI set completed...start building buffer zone\n";
    // building buffer zone 
    for (size_t r=0; r<c_hierarchy.Row; r++) {
        r0 = r*dim2;
        for (size_t c=0; c<c_hierarchy.Col; c++) {
            c0 = c*c_hierarchy.Height;
            for (size_t h=0; h<c_hierarchy.Height; h++) {
                k = r0 + c0 + h;
                l = c_hierarchy.level[k];
                if (u_map[k]==BACKGROUND) {
                    rad = R2[l-c_hierarchy.l_th];                     
                    struct cube_<int> pos_ = {r, c, h}; 
                    if(buffer_zone<T1, T2>(c_hierarchy, u_map, pos_, (int)rad)) 
                        u_map[k]=BUFFER_ZONE;
                }
            }
        }
    }
}

} // namespace mgard 


