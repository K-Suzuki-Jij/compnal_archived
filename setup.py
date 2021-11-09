from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension

ext_modules = [
    Pybind11Extension(
         "compnal",
        sorted(glob("wrapper/src/*.cpp")),  # Sort source files for reproducibility
    ),
]

setup(
    name='compnal',
    version=0.1,
    ext_modules=ext_modules
)
