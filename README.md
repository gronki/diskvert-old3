# diskvert

## Installation

### Numerical codes

Program requires **GCC 6** (or newer), including **gfortran** or other compiler suite (such as Intel Parallel Studio).
Developement package for Lapack, OpenBLAS or other compatible library is required.

```sh
cd build
make
make install
```

If one has no root access, installation can be made in user's home directory.
For example, to install in ```~/.local/bin``` and ```~/.local/lib```, instead of the last line one would run:

```sh
make install prefix=~/.local
```

Now to execute the commands freely, it's advisable to add it to system ```$PATH``` variable.
This command can be also placed in ```.bashrc``` or ```.cshrc``` file.

```sh
export PATH="$PATH:$HOME/.local/bin"
```

### Debug build

By default, program is built using compiler flags that provide the fastest execution.
If an error (such as segmentation fault) occurs, it may be helpful to build the program in a way that allows easy debugging.
This can be achieved by the ```BUILD``` switch:

```
make BUILD=debug
```

This will enable the array checking and add debug information, but will slow the execution.

### Building the library and python package

Python package, containing plotting scripts and bindings to run the code from Python has been provided.
The package needs the code in shared library format.
The shared library needs rebuilding the whole project.

```sh
make clean
make lib
# to install system-wide (requires root)
make install-lib
# to install for the current user
make install-lib prefix=~/.local
```

To install the Python package, a setup script is provided.
Is it advised (but not required) that you use a virtual environment.
Installing the package system-wide is risky and can conflict with your distribution's packages.

```sh
cd ..
virtualenv venv
. venv/bin/activate
cd python
python setup.py install
```

Alternatively, once can install the Python package for the current user only:

```sh
python setup.py install --user
```

### Other compiler vendors

Diskvert can be compiled using **icc** and **ifort**.
This is achieved by using the ```VENDOR``` switch, which will also link the program against MKL rather than standard LAPACK.

```sh
make VENDOR=intel
```

Note that PGI suite is not supported at the moment due to very poor support of the latest standard by PGI compiler.
