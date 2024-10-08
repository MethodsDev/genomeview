#!/usr/bin/env python

from setuptools import setup, find_packages, Extension

def get_version(string):
    """ Parse the version number variable __version__ from a script. """
    import re
    version_re = r"^__version__ = ['\"]([^'\"]*)['\"]"
    version_str = re.search(version_re, string, re.M).group(1)
    return version_str


setup(
    name="genomeview",
    version=get_version(open('genomeview/__init__.py').read()),
    description="genomeview",
    author="Noah Spies",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        'genomeview': ['templates/*.html'],
    },

    # setup_requires=["pypandoc"],
    install_requires=["pysam", "numpy", "pandas", "pyBigWig", "intervaltree", "Jinja2"], 
    python_requires=">=3.3",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
)
