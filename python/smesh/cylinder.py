#!/usr/bin/env python3

import gmsh
import numpy as np
import sys, getopt


def main(argv):
    usage = (
        f"usage: {argv[0]} <radius> <height> <output_db> "
        "[--refinements=N] [--order=N] "
        "[--mesh_size_min=S] [--mesh_size_max=S]"
    )

    if len(argv) < 4:
        print(usage)
        exit(1)

    radius = float(argv[1])
    height = float(argv[2])
    output = argv[3]
    mesh_size_min = None
    mesh_size_max = None
    nrefs = 0
    order = 1

    try:
        opts, args = getopt.getopt(
            argv[4:],
            "h",
            ["refinements=", "order=", "mesh_size_min=", "mesh_size_max=", "help"],
        )
    except getopt.GetoptError as err:
        print(err)
        print(usage)
        sys.exit(1)

    if args:
        print(f"unexpected positional arguments: {' '.join(args)}")
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(usage)
            sys.exit()
        elif opt == "--refinements":
            nrefs = int(arg)
        elif opt == "--order":
            order = int(arg)
        elif opt == "--mesh_size_min":
            mesh_size_min = float(arg)
        elif opt == "--mesh_size_max":
            mesh_size_max = float(arg)

    if order < 1:
        print("--order must be >= 1")
        sys.exit(1)

    gmsh.initialize(argv=["", "-bin"])
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.option.setNumber("Mesh.SaveAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm3D", 1)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
    gmsh.option.setNumber("Mesh.SecondOrderIncomplete", 0)
    gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)

    if mesh_size_min is not None:
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", mesh_size_min)

    if mesh_size_max is not None:
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", mesh_size_max)

    model = gmsh.model()
    model.add("Cylinder")
    model.setCurrent("Cylinder")
    model.occ.addCylinder(0, 0, 0, 0, height, 0, radius, tag=1)

    inlet = 1
    outlet = 2
    wall = 3

    walls = []

    for s in gmsh.model.occ.getEntities(dim=2):
        com = gmsh.model.occ.getCenterOfMass(s[0], s[1])

        if com[1] < 1e-9:
            gmsh.model.addPhysicalGroup(s[0], [s[1]], inlet)
            gmsh.model.setPhysicalName(s[1], inlet, "sinlet")
        elif com[1] >= height - 1e-9:
            gmsh.model.addPhysicalGroup(s[0], [s[1]], outlet)
            gmsh.model.setPhysicalName(s[1], outlet, "soutlet")
        else:
            walls.append(s[1])

    gmsh.model.addPhysicalGroup(2, walls, wall)

    volumes = gmsh.model.occ.getEntities(dim=3)
    model.occ.rotate(volumes, 0, 0, 0, 0, 0, 1, -np.pi / 2)
    model.occ.translate(volumes, -height / 2, 0, 0)

    # Generate mesh
    model.occ.synchronize()
    model.mesh.generate(3)

    for r in range(0, nrefs):
        model.mesh.refine()

    model.mesh.optimize("Netgen")

    if order > 1:
        model.mesh.setOrder(order)
        model.mesh.optimize("HighOrder")

    gmsh.write(output)
    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv)
