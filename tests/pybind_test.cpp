#include <iostream>
#include <pybind11/pybind11.h>

namespace py = pybind11;

void helloworld(){
    std::cout << "Hello world!" << std::endl;
    return;
}


PYBIND11_MODULE(pybind_test, m) {
    m.doc() = "pybind11 test plugin";
    m.def("helloworld", &helloworld, "The start of every program");
}

