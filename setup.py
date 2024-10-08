from importlib.metadata import entry_points
from setuptools import setup,find_packages

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext
from pybind11 import get_cmake_dir
from glob import glob
import sys

__version__ = "4.0.1"
requires_packages = [
    "numpy", "biopython", "pandas","pybind11"]
ext_modules = [
    Pybind11Extension("diverge._type2cpp",
                      sorted(glob('src/type2/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._gu99cpp",
                      sorted(glob('src/gu99/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._gu2001cpp",
                      sorted(glob('src/gu2001/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._fdrcpp",
                      sorted(glob('src/fdr/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._rvscpp",
                      sorted(glob('src/rvs/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._typeOneAnalysiscpp",
                      sorted(glob('src/typeOneAnalysis/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._asymcpp",
                      sorted(glob('src/asym/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    Pybind11Extension("diverge._effectivecpp",
                      sorted(glob('src/effective/*.c*')),
                      define_macros=[('VERSION_INFO', __version__)],
                      language='c++',
                      extra_compile_args=['-w'],
                      ),
    
]
# ext_modules[0]._add_cflags('-g')
setup(
    name="diverge",
    version=__version__,
    author="cyc",
    author_email="3170103839@zju.edu.cn",
    url="https://github.com/zjupgx/diverge4",
    description="DIVERGE python version.",
    long_description="",
    ext_modules=ext_modules,
    # extras_require={"test": "pytest"},
    packages=find_packages(),
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    install_requires = requires_packages,
)

