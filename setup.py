import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name="scripts",
    version="0.1",
    description="Install functions to run analysis to examine relationship between core and accessory gene expression",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/greenelab/core-accessory-interactome",
    author="Alexandra Lee",
    author_email="alexjlee.21@gmail.com",
    license="BSD 3-Clause",
    packages=["scripts"],
    zip_safe=False,
)
