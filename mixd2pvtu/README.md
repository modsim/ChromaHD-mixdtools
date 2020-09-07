# MIXD2PVTU

This program can convert MIXD files to PVTU for visualization. It is MPI enabled. The source code is derived from Loic Wendling's projects given out to students of ParaCompSS19 and modified to handle tetrahedrons. 

The code is relatively crude, but functional. 

# Compilation 

Prerequisites: `VTK, MPI`

```
mkdir build
cd build
cmake ..
make
```

`mpirun -np 4 settings.in`

# Notes
- [PVTU generated with MPI doesn't contain all the pieces](https://github.com/libMesh/libmesh/issues/1179)
- PVTU can be manually modified to add remaining pieces. Do a little vim-fu.
- Works only on semidiscrete meshes
- minf file nn/ne ordering matters 
- Timestepping implemented (nrec)
- Timestriding implemented (nrecstride. use ncrecstride=2 for spacetime data out) 
- Timestep lengths are not used 
- Example settings file is provided with source code

# Todos
- [ ] Error handling & safety
- [ ] Write out to subfolders
- [ ] Include simple tests
- [ ] General restructuring & improvements
- [ ] Better comments and documentation for modified files
- [ ] Fix minf file read to be order agnostic
- [ ] Allow naming scalar quantities
- [ ] Fix MPI PVTU with all pieces
- [ ] Test and ensure that data splitting happens properly. I seem to recall tiny (possibly trivial) issues in the algorithm back when I was assisting corrections.
