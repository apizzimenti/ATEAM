
from pathlib import Path
from setuptools import find_packages, setup

requirements = [
    "numpy", "galois", "matplotlib", "pandas", "scipy"
]

# Set the version --- ensure that the latest tag matches this value.
VERSION = "0.0.1"

# Description.
here = Path(__file__).parent
DESCRIPTION = (here / "README.md").read_text()

setup(
    name="potts",
    version=VERSION,
    author="Paul Duncan, Anthony E. Pizzimenti, Benjamin Schweinhart",
    author_email="apizzime@gmu.edu",
    description="Generalized Potts model.",
    long_description=DESCRIPTION,
    long_description_content_type="text/markdown",
    url="https://github.com/apizzimenti/potts",
    packages=find_packages(exclude=["test"]),
    install_requires=requirements,
    include_package_data=True,
    extras_require={
        "dev": ["pdoc3", "flake8", "pytest", "autopep8", "pytest-cov", "black", "isort"]
    },
)