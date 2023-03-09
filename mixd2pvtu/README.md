# MIXD2PVTU

This program can convert MIXD files to PVTU for visualization. It is MPI enabled. The source code is derived from Loic Wendling's projects given out to students of ParaCompSS19 and modified to handle tetrahedrons along with triangles with run-time element specification. 

The code is relatively crude, but functional. 

# Compilation 

## Using Nix

Assuming you have `nix` installed, you can just run `nix-build` in this directory to get everything installed and ready.

Check out [nix-portable](https://github.com/DavHau/nix-portable).

## Manual Build

Prerequisites: `VTK (with MPI support), MPI`

VTK can be installed with OSMESA with the following commands after downloading and unzipping the source:

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$HOME/local/VTK/9.1.0 -DVTK_USE_MPI=ON -DVTK_WRAP_PYTHON=ON -DVTK_SMP_IMPLEMENTATION_TYPE=TBB -DVTK_OPENGL_HAS_OSMESA=ON -DVTK_USE_X=OFF ..
make -j $(nprocs) install
```

Then, ensure that `$HOME/local/VTK/9.1.0/include` is in `$CPATH` and `$HOME/local/VTK/9.1.0/lib` (or `lib64`) is in `$LD_LIBRARY_PATH`.

Also, probably good to symlink `$HOME/local/VTK/9.1.0/include/vtk-9.1` -> `$HOME/local/VTK/9.1.0/include/vtk` since our includes use `vtk/` now.

Compiling MIXD2PVTU should then be simple:

```
mkdir build
cd build
cmake ../src
make install
```

Newer versions of cmake should allow the following neat one-liner:

```
cmake -B build src && make -C build install
```

# Usage
The program can be configured via a text based input file or just commandline arguments.

`mpirun -np 4 mixd2pvtu settings.in`

`mixd2pvtu -t vis -o outpath/ --ndf 2 --nrec 256`

# Notes
- Example settings file is provided with source code
- [PROB] For large meshes on single core, nnc*nsd etc is larger than INT_MAX, and overflows the int count parameters in MPI calls
- Works with VTK 9.1
- See [this](https://vtk.org/Wiki/VTK/Tutorials/CMakeListsFile) in case cmake can't find VTK.

# Todos
- [ ] Error handling & safety
- [ ] Improve allocations
- [ ] script simple tests
- [ ] General restructuring & improvements
- [ ] Better comments and documentation 
- [ ] Allow naming scalar quantities
- [ ] Test and ensure that data splitting happens properly. I seem to recall tiny (possibly trivial) issues in the mnc/nnc algorithm back when I was assisting corrections.
- [ ] Implement mixdclass to handle mixd data??
- [ ] Timing outputs for different sections of the code. Also per data write. 
- [ ] Fix MPI int count overflows for large meshes. 
- [ ] Write unit tests
- [CRIT] [ ] There should be a way to read mixd in parallel and write out non-parallel pvtu files. This is necessary for bead loading because connectivity filter doesn't work in parallel
