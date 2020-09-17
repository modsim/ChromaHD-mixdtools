# MIXD2PVTU

This program can convert MIXD files to PVTU for visualization. It is MPI enabled. The source code is derived from Loic Wendling's projects given out to students of ParaCompSS19 and modified to handle tetrahedrons. 

The code is relatively crude, but functional. 

# Compilation 

Prerequisites: `VTK (with MPI support), MPI`

```
mkdir build
cd build
cmake ..
make
```

`mpirun -np 4 settings.in`

# Notes
- [VTK 7.1 Update and VTK with MPI](https://github.com/libMesh/libmesh/issues/1179)
- Works only on semidiscrete meshes
- minf file nn/ne ordering matters 
- Timestepping implemented (nrec)
- Timestriding implemented (nrecstride. use ncrecstride=2 for spacetime data out) 
- Timestep lengths are not used 
- Example settings file is provided with source code
- ISSUE: for large meshes, it gets stuck after mxyz read.

# Todos
- [ ] Error handling & safety
- [ ] Fix memory leaks
- [X] Use long offsets and sizes
- [X] Write out to subfolders
- [ ] script simple tests
- [ ] General restructuring & improvements
- [ ] Better comments and documentation 
- [X] Clean STDOUT
- [X] Fix minf file read to be order agnostic
- [ ] Allow naming scalar quantities
- [X] Fix MPI PVTU with all pieces
- [ ] Test and ensure that data splitting happens properly. I seem to recall tiny (possibly trivial) issues in the mnc/nnc algorithm back when I was assisting corrections.
- [X] Implement make install
- [X] check that installed bin works fine
- [X] Implement filesize checks. 
- [X] Implement spacetime mesh support
- [X] Implement nrecoffset 
- [ ] Implement mixdclass to handle mixd data??
- [ ] -march, -mtune flags?
- [X] Allow stringing together multiple data files (with offset of 1 in between them)
