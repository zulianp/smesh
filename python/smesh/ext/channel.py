#!/usr/bin/env python3

import getopt
import sys
import warnings

warnings.filterwarnings("ignore")


def main(argv):
    usage = (
        f"usage: {argv[0]} <length> <width> <height> <radius> <output_db> "
        "[--cylinder_x=X] [--cylinder_z=Z] "
        "[--refinements=N] [--order=N] "
        "[--mesh_size_min=S] [--mesh_size_max=S]"
    )

    if any(arg in ("-h", "--help") for arg in argv[1:]):
        print(usage)
        sys.exit(0)

    if len(argv) < 6:
        print(usage)
        sys.exit(1)

    length = float(argv[1])
    width = float(argv[2])
    height = float(argv[3])
    radius = float(argv[4])
    output = argv[5]

    cylinder_x = length * 0.2
    cylinder_z = height * 0.5
    nrefs = 0
    order = 1
    mesh_size_min = None
    mesh_size_max = None

    try:
        opts, args = getopt.getopt(
            argv[6:],
            "h",
            [
                "cylinder_x=",
                "cylinder_z=",
                "refinements=",
                "order=",
                "mesh_size_min=",
                "mesh_size_max=",
                "help",
            ],
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
        elif opt == "--cylinder_x":
            cylinder_x = float(arg)
        elif opt == "--cylinder_z":
            cylinder_z = float(arg)
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

    import gmsh
    import numpy as np

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
    model.add("DFG 3D")
    model.setCurrent("DFG 3D")

    channel = model.occ.addBox(0, 0, 0, length, width, height)
    obstacle = model.occ.addCylinder(cylinder_x, 0, cylinder_z, 0, width, 0, radius)
    fluid = model.occ.cut([(3, channel)], [(3, obstacle)])

    model.occ.synchronize()

    volumes = model.getEntities(dim=3)
    assert volumes == fluid[0]

    fluid_marker = 11
    model.addPhysicalGroup(3, [volumes[0][1]], fluid_marker)
    model.setPhysicalName(3, fluid_marker, "Fluid volume")

    inlet_marker, outlet_marker, wall_marker, obstacle_marker = 1, 3, 5, 7
    tol = 1e-8 * max(
        1.0, length, width, height, radius, abs(cylinder_x), abs(cylinder_z)
    )
    inlet_surfaces = []
    outlet_surfaces = []
    walls = []
    obstacles = []

    boundary_surfaces = model.getBoundary(volumes, oriented=False, recursive=False)

    for dim, tag in boundary_surfaces:
        xmin, ymin, zmin, xmax, ymax, zmax = model.getBoundingBox(dim, tag)
        if np.isclose(xmin, 0.0, atol=tol) and np.isclose(xmax, 0.0, atol=tol):
            inlet_surfaces.append(tag)
        elif np.isclose(xmin, length, atol=tol) and np.isclose(xmax, length, atol=tol):
            outlet_surfaces.append(tag)
        elif (
            (np.isclose(ymin, 0.0, atol=tol) and np.isclose(ymax, 0.0, atol=tol))
            or (np.isclose(ymin, width, atol=tol) and np.isclose(ymax, width, atol=tol))
            or (np.isclose(zmin, 0.0, atol=tol) and np.isclose(zmax, 0.0, atol=tol))
            or (
                np.isclose(zmin, height, atol=tol)
                and np.isclose(zmax, height, atol=tol)
            )
        ):
            walls.append(tag)
        else:
            obstacles.append(tag)

    if not inlet_surfaces:
        raise RuntimeError("failed to identify inlet surface")

    inlet_surfaces = sorted(set(inlet_surfaces))
    outlet_surfaces = sorted(set(outlet_surfaces))
    walls = sorted(set(walls))
    obstacles = sorted(set(obstacles))

    model.addPhysicalGroup(2, inlet_surfaces, inlet_marker)
    model.setPhysicalName(2, inlet_marker, "Fluid inlet")
    if outlet_surfaces:
        model.addPhysicalGroup(2, outlet_surfaces, outlet_marker)
        model.setPhysicalName(2, outlet_marker, "Fluid outlet")
    model.addPhysicalGroup(2, walls, wall_marker)
    model.setPhysicalName(2, wall_marker, "Walls")
    model.addPhysicalGroup(2, obstacles, obstacle_marker)
    model.setPhysicalName(2, obstacle_marker, "Obstacle")

    resolution = mesh_size_min if mesh_size_min is not None else radius / 10.0
    far_resolution = mesh_size_max if mesh_size_max is not None else 20.0 * resolution
    obstacle_mid_resolution = min(far_resolution, 4.0 * resolution)
    inlet_resolution_min = min(far_resolution, 5.0 * resolution)
    inlet_mid_resolution = max(
        inlet_resolution_min, min(far_resolution, 8.0 * resolution)
    )
    inlet_resolution_max = max(
        inlet_mid_resolution, min(far_resolution, 12.0 * resolution)
    )

    distance = model.mesh.field.add("Distance")
    model.mesh.field.setNumbers(distance, "FacesList", obstacles)

    threshold = model.mesh.field.add("Threshold")
    model.mesh.field.setNumber(threshold, "IField", distance)
    model.mesh.field.setNumber(threshold, "LcMin", resolution)
    model.mesh.field.setNumber(threshold, "LcMax", obstacle_mid_resolution)
    model.mesh.field.setNumber(threshold, "DistMin", 0.5 * radius)
    model.mesh.field.setNumber(threshold, "DistMax", 1.5 * radius)

    obstacle_transition = model.mesh.field.add("Threshold")
    model.mesh.field.setNumber(obstacle_transition, "IField", distance)
    model.mesh.field.setNumber(obstacle_transition, "LcMin", obstacle_mid_resolution)
    model.mesh.field.setNumber(obstacle_transition, "LcMax", far_resolution)
    model.mesh.field.setNumber(obstacle_transition, "DistMin", 1.5 * radius)
    model.mesh.field.setNumber(obstacle_transition, "DistMax", 4.0 * radius)

    inlet_dist = model.mesh.field.add("Distance")
    model.mesh.field.setNumbers(inlet_dist, "FacesList", inlet_surfaces)

    inlet_threshold = model.mesh.field.add("Threshold")
    model.mesh.field.setNumber(inlet_threshold, "IField", inlet_dist)
    model.mesh.field.setNumber(inlet_threshold, "LcMin", inlet_resolution_min)
    model.mesh.field.setNumber(inlet_threshold, "LcMax", inlet_mid_resolution)
    model.mesh.field.setNumber(inlet_threshold, "DistMin", 0.1)
    model.mesh.field.setNumber(inlet_threshold, "DistMax", 0.35 * length)

    inlet_transition = model.mesh.field.add("Threshold")
    model.mesh.field.setNumber(inlet_transition, "IField", inlet_dist)
    model.mesh.field.setNumber(inlet_transition, "LcMin", inlet_mid_resolution)
    model.mesh.field.setNumber(inlet_transition, "LcMax", inlet_resolution_max)
    model.mesh.field.setNumber(inlet_transition, "DistMin", 0.35 * length)
    model.mesh.field.setNumber(inlet_transition, "DistMax", 0.75 * length)

    minimum = model.mesh.field.add("Min")
    model.mesh.field.setNumbers(
        minimum,
        "FieldsList",
        [threshold, obstacle_transition, inlet_threshold, inlet_transition],
    )
    model.mesh.field.setAsBackgroundMesh(minimum)

    model.mesh.generate(3)

    for _ in range(nrefs):
        model.mesh.refine()

    model.mesh.optimize("Netgen")

    if order > 1:
        model.mesh.setOrder(order)
        model.mesh.optimize("HighOrder")

    gmsh.write(output)
    gmsh.finalize()


if __name__ == "__main__":
    main(sys.argv)
