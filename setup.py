import setuptools

from modbed import version

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='modbedtools',
    version=version.__version__,
    scripts=['modbedtools'],
    author="Daofeng Li",
    author_email="lidaof@gmail.com",
    description="Generate modbed track files for visualization on WashU Epigenome Browser",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lidaof/modbedtools",
    packages=setuptools.find_packages(),
    install_requires=[
        'pysam',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
