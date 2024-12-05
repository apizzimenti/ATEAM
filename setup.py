
from setuptools import setup
from Cython.Build import cythonize
import os

#########################
### INSTALLATION NOTE ###
#########################

# PHAT is often *very difficult* to install. Before installing, please ensure
# that the following are installed and accessible in your system:
#
#   • PyBind11, setuptools, wheel (via pip)
#   • g++ (preferably >12)
#
# When installing PHAT, please use the arguments
#
#   pip install --use-deprecated=legacy-resolver --no-binary :all: phat
#
# To ensure installation.

if "CC" in os.environ: os.environ["CC"] = os.environ["CC"]

setup(
    ext_modules=cythonize(
        "ateam/**/*.pyx"
    )
)