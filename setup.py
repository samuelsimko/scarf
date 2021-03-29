import setuptools
import pybind11
from pybind11.setup_helpers import Pybind11Extension, build_ext

with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

extensions = [
    Pybind11Extension(
        "healpywrappercpp",
        sources=["./src/healpywrapper/healpywrappercpp.cc"],
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
    name="healpywrapper-samuelsimko",
    version="0.0.1",
    author="Samuel Simko",
    author_email="samuel@simko.info",
    description="A wrapper to healpy",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    url="https://github.com/samuelsimko/healpy-wrapper",
    project_urls={
        "Bug Tracker": "https://github.com/samuelsimko/healpy-wrapper/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    ext_modules=extensions,
    python_requires=">=3.6",
)
