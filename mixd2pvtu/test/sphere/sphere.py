#!/usr/bin/env python

import gmsh

gmsh.initialize()

gmsh.model.add("Sphere")

gmsh.logger.start()

radius = 5e-5

gmsh.model.occ.addSphere(0,0,0,radius,tag=1)
dt_sph = [(3,1)]

gmsh.model.occ.synchronize()

gmsh.model.mesh.setSize(gmsh.model.getEntities(0), radius * 2 * 0.1)

boundaries = gmsh.model.getBoundary(dt_sph,False,False,False)

boundary_tags = [ y for _,y in boundaries]

print("nboundaries:", len(boundary_tags))
print(boundary_tags)

gmsh.model.addPhysicalGroup(2, [boundary_tags[0]], 1)
gmsh.model.setPhysicalName(2,1,"Wall")

gmsh.model.addPhysicalGroup(3, [1], 4)
gmsh.model.setPhysicalName(3,4,"Volume")

gmsh.option.setNumber('Print.GeoOnlyPhysicals', 1)

gmsh.model.mesh.generate(3)

gmsh.write("sph.msh2")
gmsh.write("sph.vtk")
