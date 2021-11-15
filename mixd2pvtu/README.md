# MIXD2PVTU

This program can convert MIXD files to PVTU for visualization. It is MPI enabled. The source code is derived from Loic Wendling's projects given out to students of ParaCompSS19 and modified to handle tetrahedrons along with triangles with run-time element specification. 

The code is relatively crude, but functional. 

# Compilation 

Prerequisites: `VTK (with MPI support), MPI`

```
mkdir build
cd build
cmake ../src
make install
```

# Usage
The program can be configured via a text based input file or just commandline arguments.

`mpirun -np 4 mixd2pvtu settings.in`

`mixd2pvtu -t vis -o outpath/ --ndf 2 --nrec 256`

# Notes
- Example settings file is provided with source code
- [PROB] For large meshes on single core, nnc*nsd etc is larger than INT_MAX, and overflows the int count parameters in MPI calls
- Works with VTK 9.1

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
- [ ] Dynamic handling of element types.
- [ ] Rename tetMesh to just Mesh.
- [ ] Write unit tests
