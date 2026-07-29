from setuptools import setup, find_packages

setup(
    name='strand_design',
    version='0.1',
    packages=find_packages(),
    description='Given a ssOrigami,\
        design the sequence which will embbed the coding sequence in the strucutre.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/mlsample/mRNA-Design',
    author='Matthew Sample',
    author_email='matsample1@gmail.com',
    license='GPL-3.0',
    install_requires=[
    'numpy',
    'pandas',
    'matplotlib',
    'seaborn',
    'biopython',
    'pytest'
    # ipy_oxdna, oxDNA_analysis_tools/oxpy and nupack are not on PyPI
    # and must be installed separately
    ],
)
