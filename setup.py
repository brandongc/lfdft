# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()


setup(
    name='lfdft',
    version='0.1',
    description='Density functional theory with lagrange function basis',
    long_description=long_description,
    url='https://github.com/brandongc/lfdft',
    author='Brandon Cook',
    author_email='brandongc@gmail.com',
    license='Apache Software License',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 2.7'
    ],
    keywords='dft lagrange',
    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
    install_requires=['numpy', 'scipy'],
    package_data={
        'lfdft' : ['data/*.dat']
    },
    scripts=['bin/lfdft'],
)
