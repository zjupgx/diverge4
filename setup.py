import sys
from glob import glob

from setuptools import setup, find_packages
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11 import get_cmake_dir

__version__ = "4.0.1"

with open("requirements.txt", "r") as f:
    REQUIRED_PACKAGES = f.read().splitlines()

def get_extension_modules():
    """Define and return the list of extension modules."""
    common_compile_args = ['-w']
    common_macros = [('VERSION_INFO', __version__)]

    extension_modules = [
        Pybind11Extension(
            f"diverge._{module}cpp",
            sorted(glob(f'src/{module}/*.c*')),
            define_macros=common_macros,
            language='c++',
            extra_compile_args=common_compile_args,
        )
        for module in [
            "type2", "gu99", "gu2001", "fdr", "rvs",
            "typeOneAnalysis", "asym", "effective"
        ]
    ]

    return extension_modules

setup(
    name="diverge",
    version=__version__,
    author="cyc",
    author_email="3170103839@zju.edu.cn",
    url="https://github.com/zjupgx/diverge4",
    description="DIVERGE 4.0: A Python package for detecting functional divergence in protein family evolution",
    long_description="""DIVERGE 4.0 is a comprehensive software package for detecting functional divergence in protein family evolution. 
    It implements various statistical methods to identify amino acid sites that have undergone functional divergence, 
    including Type-I and Type-II functional divergence. This package is an essential tool for researchers 
    in evolutionary biology, bioinformatics, and computational biology.""",
    long_description_content_type="text/plain",
    ext_modules=get_extension_modules(),
    packages=find_packages(),
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    install_requires=REQUIRED_PACKAGES,
    keywords="bioinformatics evolution functional-divergence protein-family",
)

