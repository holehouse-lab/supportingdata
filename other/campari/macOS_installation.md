# CAMPARI 5.0 Installation Guide for macOS (Apple Silicon)

#### Created 2026-02-26

This guide documents the complete installation process for CAMPARI 5.0 on macOS with Apple Silicon (M3/M2/M1). It includes all necessary patches and workarounds for compatibility issues.

## Prerequisites

### 1. Install Xcode Command Line Tools

```bash
xcode-select --install
```

### 2. Install Homebrew (if not already installed)

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### 3. Install GFortran Compiler

```bash
brew install gcc
```

This installs gfortran as part of the GCC package. Verify the installation:

```bash
gfortran --version
```

### 4. Install Required Libraries

```bash
brew install fftw netcdf netcdf-fortran
```

**Note:** We use the macOS Accelerate framework for LAPACK/BLAS instead of OpenBLAS to avoid OpenMP ABI conflicts between LLVM and GNU OpenMP implementations.

## Source Code Modifications

Before compilation, several source files need to be patched for macOS compatibility.

### 1. Fix XDR Library for macOS

The XDR library uses Linux-specific headers and paths that don't exist on macOS.

#### a) Remove `malloc.h` include (macOS uses `stdlib.h` instead)

Edit `source/xdr/libxdrf.m4` and `source/xdr/libxdrf.c`:

Remove the line:
```c
#include <malloc.h>
```

#### b) Add `string.h` header

Add after the other includes:
```c
#include <string.h>
```

#### c) Fix include order for `xdrstdio_create`

On macOS, `stdio.h` must be included BEFORE `rpc/xdr.h` for `xdrstdio_create` to be declared. Change the include order to:

```c
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <rpc/rpc.h>
#include <rpc/xdr.h>
#include "xdrf.h"
```

#### d) Fix `ctest.c` similarly

Edit `source/xdr/ctest.c` and remove the `#include <malloc.h>` line.

#### e) Fix path to `true` command in XDR Makefile

Edit `source/xdr/Makefile` and change:
```makefile
RMCMD	= /bin/true
```
to:
```makefile
RMCMD	= true
```

(On macOS, `true` is at `/usr/bin/true`, not `/bin/true`)

## Configuration

### 1. Navigate to source directory

```bash
cd campari/source
```

### 2. Create Makefile.local

Create a file called `Makefile.local` in the source directory with the following content:

```makefile
# Makefile.local for CAMPARI on macOS (Apple Silicon/M3)
# Full installation with FFTW, NetCDF, and LAPACK support

# Example: CAMPARI_HOME=/Users/alex/tools/campari5/campari_v5_04192024/campari
CAMPARI_HOME=/path/to/your/campari

# Architecture identifier
ARCH=arm64

# Directory layout
BIN_DIR=${CAMPARI_HOME}/bin
LIB_DIR=${CAMPARI_HOME}/lib
SRC_DIR=${CAMPARI_HOME}/source

# Fortran compiler (GNU)
FF=gfortran
LFLAGS=-O3 -fopenmp

# Override default compiler settings for macOS/gfortran compatibility
# Using -std=gnu instead of -std=f2008 to avoid DATE_AND_TIME issues with newer gfortran
COMPDEFAULTS=-cpp -O3 -fall-intrinsics -ffpe-trap=invalid,zero -fdefault-real-8 -fdefault-double-8 -std=gnu

# Include paths for external libraries
NETCDF_INC=/opt/homebrew/opt/netcdf-fortran/include
FFTW_INC=/opt/homebrew/include

# Library paths
NETCDF_LIB=/opt/homebrew/opt/netcdf-fortran/lib/libnetcdff.a /opt/homebrew/lib/libnetcdf.a
FFTW_LIB=/opt/homebrew/lib/libfftw3.a /opt/homebrew/lib/libfftw3_threads.a
# Use macOS Accelerate framework for LAPACK/BLAS (avoids OpenMP ABI conflicts with OpenBLAS)
LAPACK_LIB=-framework Accelerate

# Linker: Additional libraries  
# Note: -lz, -lhdf5, and -lxml2 needed for NetCDF on macOS
EXTRA_LIBS=${FFTW_LIB} ${LAPACK_LIB} ${NETCDF_LIB} -L/opt/homebrew/lib -lhdf5 -lhdf5_hl -lz -lcurl -lxml2

# Compiler flags 
# -DLINK_NETCDF : Enable NetCDF trajectory I/O
# -DLINK_XDR    : Enable XTC trajectory I/O (compressed GROMACS format)
# -DLINK_LAPACK : Enable LAPACK for constraint solving and PCA
# -DLINK_FFTW   : Enable FFTW for particle-mesh Ewald
# -DDISABLE_ERFTAB : Use compiler intrinsics for erf/erfc
# -DDISABLE_FLOAT  : Disable single precision (required)
FFLAGS=-DLINK_NETCDF -DLINK_XDR -DLINK_LAPACK -DDISABLE_ERFTAB -DDISABLE_FLOAT -DLINK_FFTW -I${FFTW_INC} -I${NETCDF_INC}

# MPI settings (if using MPI - optional)
MPIFF=mpif90
MPILFLAGS=-O3
MPIFFLAGS=${FFLAGS} 
MPIEXTRA_LIBS=${EXTRA_LIBS}

# OpenMP threading support
THREADFLAGS=-fopenmp

# C compiler for XDR library
CC=gcc
XDRFFLAGS=-O3
XDRCFLAGS=-O3

# C++ compiler (for PostgreSQL module if needed)
CXX=g++
PSQLCXXFLAGS=-O3
PSQLEXTRALIBS=

# HSL library settings (disabled by default - requires separate license)
HSLFFLAGS=${XDRFFLAGS}
HSLEXTRALIBS=-framework Accelerate
```

**Important:** Replace `/path/to/your/campari` with the actual path to your CAMPARI installation.

### 3. Create build directories

```bash
mkdir -p ../lib/arm64/threads ../bin/arm64
```

### 4. Generate dependencies file

```bash
./make_dependencies.sh .
```

## Compilation

### Build the serial version

```bash
make -f Makefile.manual campari
```

Compilation will take several minutes. You'll see warnings about deprecated C function prototypes in the XDR library - these can be safely ignored.

### Verify the build

```bash
ls -la ../bin/arm64/campari
```

You should see an executable file of approximately 9 MB.

### Test the executable

```bash
../bin/arm64/campari --help
```

You should see the CAMPARI welcome message and version information.

## Optional: Build Thread-Parallel Version

To build the OpenMP-enabled version:

```bash
make -f Makefile.manual campari_threads
```

The executable will be created at `bin/arm64/campari_threads`.

## Post-Installation

### Add to PATH

Add the following to your `~/.zshrc`:

```bash
export CAMPARI_HOME="/Users/alex/tools/campari5/campari_v5_04192024/campari"
export PATH="$CAMPARI_HOME/bin/arm64:$PATH"
```

Then reload:

```bash
source ~/.zshrc
```

### Verify installation

```bash
campari --help
```

## Features Enabled

With this configuration, the following features are enabled:

| Feature | Description |
|---------|-------------|
| FFTW | Fast Fourier Transforms for particle-mesh Ewald summation |
| NetCDF | Reading/writing NetCDF trajectory files |
| LAPACK | Linear algebra for constraint solving and PCA analysis |
| XDR | Reading/writing XTC trajectory files (GROMACS format) |

## Features Not Enabled

The following optional features are not enabled by default:

| Feature | Reason |
|---------|--------|
| HSL | Requires separate license from STFC |
| PostgreSQL | Requires libpqxx and PostgreSQL development libraries |
| Python interface | Requires additional configuration |
| MPI | Requires MPI installation (OpenMPI or MPICH) |

## Troubleshooting

### Error: `'malloc.h' file not found`

This means the XDR library patches weren't applied. See the "Fix XDR Library for macOS" section above.

### Error: `xdrstdio_create undeclared`

The include order in the XDR library is wrong. Make sure `stdio.h` is included before `rpc/xdr.h`.

### Error: `__kmpc_*` symbols not found

This occurs when mixing LLVM OpenMP (from OpenBLAS) with GNU OpenMP (from gfortran). Use the macOS Accelerate framework instead of OpenBLAS:

```makefile
LAPACK_LIB=-framework Accelerate
```

### Error: `DATE_AND_TIME` non-default kind

This occurs with `-std=f2008` flag in newer gfortran versions. Use `-std=gnu` instead in COMPDEFAULTS.

### Error: `/bin/true: No such file or directory`

macOS has `true` at `/usr/bin/true`. Edit `source/xdr/Makefile` and change `/bin/true` to `true`.

## References

- [CAMPARI Official Documentation](https://campari.sourceforge.net/V5/documentation.html)
- [CAMPARI Installation Guide](https://campari.sourceforge.net/V5/install.html)
- [FFTW Library](http://www.fftw.org/)
- [NetCDF Library](https://www.unidata.ucar.edu/software/netcdf/)
