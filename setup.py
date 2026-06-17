# BioMark - Package Setup
# Author: Amir Shahbazi

from setuptools import setup, find_packages
import os

# Read long description from README
here = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(here, "README.md"), encoding="utf-8") as f:
      long_description = f.read()

setup(
      name="biomark",
      version="0.5.0",
      author="Amir Shahbazi",
      author_email="",
      description="Open-source bioinformatics benchmark tool for hardware evaluation",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url="https://github.com/shahbazigenomics/BioMark",
      project_urls={
                "Bug Tracker": "https://github.com/shahbazigenomics/BioMark/issues",
                "Documentation": "https://github.com/shahbazigenomics/BioMark#readme",
      },
      license="MIT",
      classifiers=[
                "Development Status :: 4 - Beta",
                "Intended Audience :: Science/Research",
                "Topic :: Scientific/Engineering :: Bio-Informatics",
                "License :: OSI Approved :: MIT License",
                "Programming Language :: Python :: 3",
                "Programming Language :: Python :: 3.8",
                "Programming Language :: Python :: 3.9",
                "Programming Language :: Python :: 3.10",
                "Programming Language :: Python :: 3.11",
                "Programming Language :: Python :: 3.12",
                "Operating System :: MacOS",
                "Operating System :: POSIX :: Linux",
      ],
      keywords=[
                "bioinformatics", "benchmark", "genomics", "hardware",
                "WES", "RNA-seq", "scRNA-seq", "metagenomics", "assembly"
      ],
      packages=find_packages(where="src"),
      package_dir={"": "src"},
      python_requires=">=3.8",
      install_requires=[
                "psutil>=5.9.0",
                "py-cpuinfo>=9.0.0",
                "numpy>=1.24.0",
                "scipy>=1.10.0",
      ],
      extras_require={
                "dev": [
                              "pytest>=7.4.0",
                              "pytest-cov>=4.1.0",
                ],
      },
      entry_points={
                "console_scripts": [
                              "biomark=main:main",
                ],
      },
      include_package_data=True,
)
