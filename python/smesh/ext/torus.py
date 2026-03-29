#!/usr/bin/env python3

import gmsh
import meshio
import numpy as np
import sys, getopt


def main(argv):
    usage = (
        f"usage: {argv[0]} <major_radius> <minor_radius> <output_db> "
        "[--refinements=N] [--order=N] "
        "[--mesh_size_min=S] [--mesh_size_max=S]"
    )

    if len(argv) < 4:
        print(usage)
        exit(1)

    major_radius = float(argv[1])  # Major radius (distance from center to tube center)
    minor_radius = float(argv[2])  # Minor radius (tube radius)
    output = argv[3]
    nrefs = 0  # Number of refinements
    order = 1
    mesh_size_min = None
    mesh_size_max = None

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
    gmsh.option.setNumber("General.Terminal", 1)
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
    model.add("Torus")
    model.setCurrent("Torus")

    # Create torus using OpenCASCADE kernel
    # addTorus(x, y, z, r1, r2, angle, tag=-1)
    # r1 = minor radius, r2 = major radius
    torus = model.occ.addTorus(0, 0, 0, major_radius, minor_radius)

    # Synchronize before meshing
    model.occ.synchronize()

    # Generate 3D mesh
    model.mesh.generate(3)

    for r in range(0, nrefs):
        print(f"refinement {r} ...")
        model.mesh.refine()

    model.mesh.optimize("Netgen")

    if order > 1:
        model.mesh.setOrder(order)
        model.mesh.optimize("HighOrder")

    # Write mesh
    gmsh.write(output)

    # Print mesh statistics
    print(f"\nMesh statistics:")
    print(f"  Major radius: {major_radius}")
    print(f"  Minor radius: {minor_radius}")
    print(f"  Refinements: {nrefs}")
    print(f"  Order: {order}")

    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv)
