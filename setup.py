from setuptools import setup, find_packages
import glob
import os
import pkg_resources

__version__ = "0.1.0"

setup(name='virataxm',
      version=__version__,
      packages=find_packages(),
      package_data={"virataxm":["ref.zip", "test/*"]},
      description='ViraTaxM: a Viral Taxonomic classification pipeline for Metagenomic sequences',
      keywords=['Bioconda Virus TaxonomicClassification Metagenomics Genus'],
      classifiers=[],
      url='https://github.com/GreyGuoweiChen/ViraTaxM.git',
      author='CHEN Guowei',
      author_email='gwchen3-c@my.cityu.edu.hk',
      license='MIT',
      install_requires=[
        'numpy>=1.23.5',
        'pandas>=2.0.3',
        'scikit-learn>=1.3.2',
        'biopython>=1.83',
    ],
      extras_require={
        "conda_only": ["mcl", "prodigal", "diamond"]
      },
      entry_points={
        "console_scripts":[
        "virataxm = virataxm.main:main"
        ]},
      include_package_data=True,
      zip_safe=False
)
