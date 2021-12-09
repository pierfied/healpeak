//
// Created by pierfied on 12/9/21.
//

#include "healpeak.h"
#include <healpix_cxx/healpix_base.h>

py::array_t<long> peak_counts(py::array_t<double, py::array::c_style | py::array::forcecast> maps,
                              py::array_t<double, py::array::c_style | py::array::forcecast> bins) {
    if (maps.ndim() == 1)
        maps.resize({1, (int) maps.shape(0)});

    int nmaps = maps.shape(0);
    int npix = maps.shape(1);
    int nbins = bins.shape(0) - 1;

    auto peak_counts = py::array_t<long>({nmaps, nbins});

    int nside = Healpix_Base::npix2nside(npix);
    auto map_base = Healpix_Base(nside, NEST, SET_NSIDE);

    double *map_data = (double *) maps.request().ptr;
    double *bin_data = (double *) bins.request().ptr;
    long *peak_count_data = (long *) peak_counts.request().ptr;
    for (int i = 0; i < nmaps * nbins; ++i) {
        peak_count_data[i] = 0;
    }

    for (int i = 0; i < npix; ++i) {
        auto neighbors = fix_arr<int, 8>();
        map_base.neighbors(i, neighbors);

        for (int j = 0; j < nmaps; ++j) {
            int ind_pix = i + j * npix;

            bool is_peak = true;
            for (int k = 0; k < 8; ++k) {
                if (neighbors[k] < 0)
                    continue;

                int ind_neighbor = neighbors[k] + j * npix;
                if (map_data[ind_pix] < map_data[ind_neighbor]) {
                    is_peak = false;
                    break;
                }
            }

            if (is_peak) {
                for (int k = 0; k < nbins; ++k) {
                    if (map_data[ind_pix] >= bin_data[k] && map_data[ind_pix] < bin_data[k + 1]) {
                        peak_count_data[k + j * nbins] += 1;
                        break;
                    }
                }
            }
        }
    }

    return peak_counts;
}

py::array_t<long> void_counts(py::array_t<double, py::array::c_style | py::array::forcecast> maps,
                              py::array_t<double, py::array::c_style | py::array::forcecast> bins) {
    if (maps.ndim() == 1)
        maps.resize({1, (int) maps.shape(0)});

    int nmaps = maps.shape(0);
    int npix = maps.shape(1);
    int nbins = bins.shape(0) - 1;

    auto void_counts = py::array_t<long>({nmaps, nbins});

    int nside = Healpix_Base::npix2nside(npix);
    auto map_base = Healpix_Base(nside, NEST, SET_NSIDE);

    double *map_data = (double *) maps.request().ptr;
    double *bin_data = (double *) bins.request().ptr;
    long *void_count_data = (long *) void_counts.request().ptr;
    for (int i = 0; i < nmaps * nbins; ++i) {
        void_count_data[i] = 0;
    }

    for (int i = 0; i < npix; ++i) {
        auto neighbors = fix_arr<int, 8>();
        map_base.neighbors(i, neighbors);

        for (int j = 0; j < nmaps; ++j) {
            int ind_pix = i + j * npix;

            bool is_void = true;
            for (int k = 0; k < 8; ++k) {
                if (neighbors[k] < 0)
                    continue;

                int ind_neighbor = neighbors[k] + j * npix;
                if (map_data[ind_pix] > map_data[ind_neighbor]) {
                    is_void = false;
                    break;
                }
            }

            if (is_void) {
                for (int k = 0; k < nbins; ++k) {
                    if (map_data[ind_pix] >= bin_data[k] && map_data[ind_pix] < bin_data[k + 1]) {
                        void_count_data[k + j * nbins] += 1;
                        break;
                    }
                }
            }
        }
    }

    return void_counts;
}
