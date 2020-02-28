from setuptools import setup, find_packages
setup(
    name='read-ICESat-2',
    version='1.0.0.6',
    description='Tools for reading and using data from the NASA ICESat-2 mission',
    url='https://github.com/tsutterley/read-ICESat-2',
    author='Tyler Sutterley',
    author_email='tsutterl@uw.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Physics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='ICESat-2 laser altimetry, ATLAS',
    packages=find_packages(),
    install_requires=['numpy','scipy','h5py','future','lxml','paramiko','scp'],
)
