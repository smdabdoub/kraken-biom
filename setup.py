# -*- coding: utf-8 -*-
from __future__ import with_statement
from setuptools import setup
 
 
def get_version():
    with open('kraken_biom.py') as f:
        for line in f:
            if line.startswith('__version__'):
                return eval(line.split('=')[-1])
 
 
def get_long_description():
    descr = []
    for fname in 'README.rst', 'CHANGELOG.rst':
        with open(fname) as f:
            descr.append(f.read())
    return '\n\n'.join(descr)
 
install_requires = ["biom-format >= 2.1.5"]

 
setup(
    name='kraken-biom',
    version=get_version(),
    description="Create BIOM-format tables from Kraken output.",
    long_description=get_long_description(),
    keywords='kraken, BIOM, metagenomics, bioinformatics, taxonomy',
    author='Shareef M. Dabdoub',
    author_email='dabdoub.2@osu.edu',
    url='https://github.com/smdabdoub/kraken-biom',
    license='MIT',
    py_modules=['kraken_biom'],
    namespace_packages=[],
    include_package_data=True,
    zip_safe=False,
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'kraken-biom = kraken_biom:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
)