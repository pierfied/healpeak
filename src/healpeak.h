//
// Created by pierfied on 12/9/21.
//

#ifndef HEALPEAK_HEALPEAK_H
#define HEALPEAK_HEALPEAK_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

py::array_t<long> peak_counts(py::array_t<double, py::array::c_style | py::array::forcecast> maps,
                              py::array_t<double, py::array::c_style | py::array::forcecast> bins,
                              std::optional<py::array_t<bool, py::array::c_style | py::array::forcecast>> mask);

py::array_t<long> void_counts(py::array_t<double, py::array::c_style | py::array::forcecast> maps,
                              py::array_t<double, py::array::c_style | py::array::forcecast> bins,
                              std::optional<py::array_t<bool, py::array::c_style | py::array::forcecast>> mask);

PYBIND11_MODULE(healpeak, m) {
    m.doc() = "Fast peak/void counting for HEALPix maps.";
    m.def("peak_counts", &peak_counts, "Compute peak counts.", py::arg("maps"), py::arg("bins"),
          py::arg("mask") = py::none());
    m.def("void_counts", &void_counts, "Compute void counts.", py::arg("maps"), py::arg("bins"),
          py::arg("mask") = py::none());
}

#endif //HEALPEAK_HEALPEAK_H
