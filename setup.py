from setuptools import setup, find_packages

setup(
    name='strand_design',
    version='0.1',
    packages=find_packages(),
    description='Given a ssOrigami,\
        design the sequence which will embbed the coding sequence in the strucutre.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/mlsample/nanostructred_expressing_nucleic_acids',
    author='Matthew Sample',
    author_email='matsample1@gmail.com',
    license='MIT',
    install_requires=[
    'numpy',
    'pytest'
    # ... other dependencies ...
    ],
    # dependencies can be listed under install_requires
)
