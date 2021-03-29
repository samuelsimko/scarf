import setuptools
import pybind11
from pybind11.setup_helpers import Pybind11Extension

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

extensions = [
    Pybind11Extension(
        "scarfcpp",
        sources=["./scarf/scarfcpp.cc"],
        include_dirs=["./ducc/src/", pybind11.get_include(True), pybind11.get_include(False)],
        extra_compile_args=[
            "-std=c++17",
            "-O3",
            "-g0",
            "-shared",
            "-march=native",
            "-Wall",
        ],
        language="c++",
    )
]

setuptools.setup(
    name="scarf",
    version="0.0.1",
    author="Samuel Simko",
    author_email="samuel@simko.info",
    description="Spherical Harmonic Transforms for CMB lensing",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/samuelsimko/scarf",
    project_urls={
        "Bug Tracker": "https://github.com/samuelsimko/scarf/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(),
    ext_modules=extensions,
    python_requires=">=3.6",
)
