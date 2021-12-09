//
// Created by pierfied on 12/9/21.
//

#ifndef HEALPEAK_HEALPEAK_H
#define HEALPEAK_HEALPEAK_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<long> peak_counts(py::array_t<double, py::array::c_style | py::array::forcecast> maps,
                              py::array_t<double, py::array::c_style | py::array::forcecast> bins);

py::array_t<long> void_counts(py::array_t<double, py::array::c_style | py::array::forcecast> maps,
                              py::array_t<double, py::array::c_style | py::array::forcecast> bins);

PYBIND11_MODULE(healpeak, m) {
    m.doc() = "Healpix peak counts";
    m.def("peak_counts", &peak_counts, "Healpix peak counts");
    m.def("void_counts", &void_counts, "Healpix void counts");
}

#endif //HEALPEAK_HEALPEAK_H
