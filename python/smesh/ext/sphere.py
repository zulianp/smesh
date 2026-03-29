#!/usr/bin/env python3

import gmsh
import meshio
import numpy as np
import sys, getopt


def main(argv):
    usage = (
        f"usage: {argv[0]} <radius> <output_db.vtk> "
        "[--refinements=N] [--order=N] "
        "[--mesh_size_min=S] [--mesh_size_max=S]"
    )

    if len(argv) < 3:
        print(usage)
        exit(1)

    radius = float(argv[1])
    output = argv[2]
    nrefs = 1
    order = 1
    mesh_size_min = None
    mesh_size_max = None

    try:
        opts, args = getopt.getopt(
            argv[3:],
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
        elif opt in ("--refinements"):
            nrefs = int(arg)
        elif opt in ("--order"):
            order = int(arg)
        elif opt in ("--mesh_size_min"):
            mesh_size_min = float(arg)
        elif opt in ("--mesh_size_max"):
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
    model.add("Sphere")
    model.setCurrent("Sphere")
    model.occ.addSphere(0, 0, 0, radius)

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
