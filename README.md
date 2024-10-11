
# What is Loci
Loci (pronounced low-sigh) is a both a C++ library and a programming framework
specifically designed for developing computational simulations of physical fields,
such as computational fluid dynamics. One advantage the framework provides is
automatic parallelization. Once an application is described within the Loci
framework, the application can be executed in parallel without change, even
though the description within the framework included no explicit parallel
directives. A particular advantage of the programming framework is that it
provides a formal framework for the development of simulation knowledge-bases
using *logic-relational* abstractions. While the approach will probably be
alien to most who begin to use it, the programming model is extremely powerful
and worth the patience required to adjust to a new way of thinking about programming.

# Compiling Loci
## Requirements
* OpenMPI
* HDF5
* GKlib (submodule)
* METIS (submodule)
* ParMETIS (submodule)
* PETSc (optional)
* LAPACK (optional)

## Quickstart
Once OpenMPI and HDF5 are compiled and added to the terminal environment (usually
accomplished through environment module files) navigate to the directory where
the Loci repo will be located and clone with recursive submodules. Next, execute
the `configure` script to setup the build directory and files. `configure --help`
will provide the most common options. If no build directory is provided as an
argument the `OBJ` directory will automatically be created to store the files
created during compilation. The `OBJ` directory will be referenced as
`$LOCI_SRC` in the Makefiles. Below is a set of example commands to clone,
configure, and install Loci.

```bash
cd $HOME/code
git clone --recurse-submodules https://github.com/mhunt0/loci.git
cd loci
./configure --prefix=$HOME/local/loci --with-mpi=${MPI_DIRECTORY}
make install
```

The result of the above commands will be installed in the user-defined *prefix*
destination directory. This directory should be referenced as `$LOCI_BASE`
when linking and compiling other software built on this Loci framework.

## Configuration Files
The `${LOCI_BASE}` directory will contain three configuration files:
* `Loci.conf` - Used to create the GNU Make compiler flags based on the `comp.conf`
  and `sys.conf`.
* `comp.conf` - (Compiler config) Used to create compiler specific flags
  (i.e. gcc, icc, etc.) that are compatible with Loci during the build process.
* `sys.conf` - (System config) Create by the `configure` script that will
  automatically find the needed libraries and executables from the user's
  `$PATH` and `$LD_LIBRARY_PATH` environment variables.


## For Developers
When working on different components/models within Loci, all of the different
libraries and executables can be compiled from the root directory. Below is a
table that describes the theme for the `Makefiles` that are consistent for each
library/executable:

| Command                       | Description                                                  |
| -----------                   | -----------                                                  |
| `make <target>`               | Compile the target                                           |
| `make install_<target>`       | Copy files to `$LOCI_BASE`                                   |
| `make clean_install_<target>` | Remove files from `$LOCI_BASE`                               |
| `make clean_<target>`         | Remove object, binary, and shared lib files from `$LOCI_SRC` |
| `make distclean_<target>`     | Remove object, binary, shared lib, and dependency files from `$LOCI_SRC` |
| `make spotless_<target>`      | Remove everything from `$LOCI_SRC` and `$LOCI_BASE`          |
| `make devhelp`                | Display environment/compilation variables                    |



# Disclaimer

The source code in the sprng is the parallel random number generator library
developed at the Florida State University. It is provided here as a convenience.
The library can be found at http://sprng.cs.fsu.edu/

The source code under FVMtools/libadf is the ADF library is an open
source library that is used by some grid convertes found in FVMtools
and is provided here as a convenience.

All other source code in this directory is provided by Mississippi
State University as free software under the Lesser GNU General Public
License version 3.  See the file [](COPYING.LESSER) for license terms.

Ed Luke  (luke@cse.msstate.edu) | Professor of Computer Science and Engineering |
Mississippi State University
