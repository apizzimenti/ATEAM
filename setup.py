
from setuptools import setup
from Cython.Build import cythonize
import os

os.environ["CC"] = "gcc-14"

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
        "potts/**/*.pyx"
    )
)