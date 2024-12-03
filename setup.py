
from setuptools import setup
from Cython.Build import cythonize
import os


if "CC" in os.environ:
    os.environ["CC"] = os.environ["CC"]

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


setup(
    # name="potts",
    # version=VERSION,
    # author="Paul Duncan, Anthony E. Pizzimenti, Benjamin Schweinhart",
    # author_email="apizzime@gmu.edu",
    # description="Generalized Potts model.",
    # long_description=DESCRIPTION,
    # long_description_content_type="text/markdown",
    # url="https://github.com/apizzimenti/potts",
    # packages=find_packages(exclude=["test"]),
    # install_requires=requirements,
    # include_package_data=True,
    # extras_require={
    #     "dev": ["pdoc3", "flake8", "pytest", "autopep8", "pytest-cov", "black", "isort"]
    # },
    ext_modules=cythonize(
        "ateam/**/*.pyx"
    )
)