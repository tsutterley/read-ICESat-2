Installation
============

Presently read-ICESat-2 is only available for use as a [GitHub repository](https://github.com/tsutterley/read-ICESat-2).
The contents of the repository can be download as a [zipped file](https://github.com/tsutterley/read-ICESat-2/archive/main.zip)  or cloned.
To use this repository, please fork into your own account and then clone onto your system.  
```bash
git clone https://github.com/tsutterley/read-ICESat-2.git
```
Can then install using `setuptools`
```bash
python setup.py install
```
or `pip`
```bash
python3 -m pip install --user .
```
Alternatively can install the utilities directly from GitHub with `pip`:
```
python3 -m pip install --user git+https://github.com/tsutterley/read-ICESat-2.git
```
Executable versions of this repository can also be tested using [Binder](https://mybinder.org/v2/gh/tsutterley/read-ICESat-2/main) and [Pangeo](https://binder.pangeo.io/v2/gh/tsutterley/read-ICESat-2/main).
