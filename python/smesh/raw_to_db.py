#!/usr/bin/env python3

import meshio
import numpy as np
import sys, getopt
import os
import glob
import pdb
from common.utils import detect_files, extension_to_dtype, extension
from common.hex27_ordering import (
    exodus_hex27_to_vtk_hex27,
    proteus_hex27_to_exodus_hex27 as proteus_hex27_to_hexahedron27,
)

import inspect

# Get the current frame
frame = inspect.currentframe()

tet4_names = ("tetra", "tetrahedron", "tetra4", "tet4", "TET4", "TETRA")
tet10_names = ("tetra10", "TET10")
hex8_names = ("hexahedron", "hex", "hex8", "HEX8", "HEX")
hex27_names = ("hexahedron27", "hex27", "HEX27")
proteus_hex27_names = ("proteus_hex27", "PROTEUS_HEX27")

quad4_names = ("quad", "quad4", "QUAD4", "QUAD")
quad9_names = ("quad9", "QUAD9", "quadshell9", "QUADSHELL9")
tri3_names = ("tri", "tri3", "TRI3", "TRI")
tri6_names = ("triangle6", "TRI6")

try:
    geom_t
except NameError:
    print("raw_to_db: self contained mode")
    geom_t = np.float32
    idx_t = np.int32

max_nodes_x_element = 27

EXODUS_OUTPUT_EXTENSIONS = (".e", ".exo", ".ex2")


def is_exodus_output(path):
    return os.path.splitext(path)[1].lower() in EXODUS_OUTPUT_EXTENSIONS


def read_simple_meta(path):
    meta = {}
    if not os.path.exists(path):
        return meta

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("- "):
                continue
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            meta[key.strip()] = value.strip()
    return meta


def element_type_from_meta(folder):
    meta_path = os.path.join(folder, "meta.yaml")
    meta = read_simple_meta(meta_path)
    for key in ("element_type", "cell_type", "elem_type"):
        if key in meta and meta[key]:
            return meta[key], meta, meta_path
    if "elem_num_nodes" not in meta:
        return None, meta, meta_path
    return None, meta, meta_path


class Field:
    def __init__(self, path):
        self.path = path
        self.type = extension_to_dtype(extension(path))
        self.data = np.fromfile(path, dtype=self.type)

        if self.type == np.int32 or self.type == np.uint8:
            print(
                f"Warning converting field from {self.type} to float32 to work with Paraview"
            )
            self.data = self.data.astype(np.float32)
            self.type = np.float32

        self.name = os.path.splitext(os.path.basename(path))[0]

    def __len__(self):
        return len(self.data)

    def __getitem__(self, index):
        return self.data[index]

    def __setitem__(self, index, value):
        self.data[index] = value

    def __str__(self):
        return f"Field(path={self.path}, type={self.type}, data={self.data})"

    def __repr__(self):
        return self.__str__()

    def check_len(self, check_len):
        if len(self.data) != check_len:
            frame = inspect.currentframe()
            print(
                f"Error in {__file__} at line {frame.f_lineno}:\n .... data length is different from number of nodes {len(self.data)} != {check_len}"
            )

            sys.exit(1)

        print(
            f"field: {self.name}, min={np.min(self.data)}, max={np.max(self.data)}, sum={np.sum(self.data)} type={self.type}"
        )


def write_transient_data(
    output_path,
    points,
    cells,
    point_data,
    cell_data,
    n_time_steps,
    time_whole,
    time_step_format,
):

    with meshio.xdmf.TimeSeriesWriter(output_path) as writer:
        writer.write_points_cells(points, cells)
        cell_data_steps = [None] * n_time_steps

        if cell_data:
            paths = cell_data.split(",")

            for p in paths:
                cell_data_files = glob.glob(p, recursive=False)
                cell_data_files.sort()

                if len(cell_data_files) < n_time_steps:
                    print(
                        f"Invalid sequence length {len(cell_data_files)} for pattern {p}"
                    )

                for i in range(0, n_time_steps):
                    if not cell_data_steps[i]:
                        cell_data_steps[i] = []

                    cell_data_steps[i].append(cell_data_files[i])

        point_data_steps = [None] * n_time_steps

        if point_data:
            paths = point_data.split(",")

            for p in paths:
                point_data_files = glob.glob(p, recursive=False)
                point_data_files.sort()

                if len(point_data_files) < n_time_steps:
                    print(
                        f"Invalid sequence length {len(point_data_files)} for pattern {p}"
                    )

                for i in range(0, n_time_steps):
                    if not point_data_steps[i]:
                        point_data_steps[i] = []

                    if i >= len(point_data_files):
                        print(
                            f"Trying to acess file that does not exists at index {i}/{n_time_steps}"
                        )

                    point_data_steps[i].append(point_data_files[i])

        for t in range(0, n_time_steps):
            name_to_point_data = {}
            name_to_cell_data = {}

            cds = point_data_steps[t]

            has_point_data = False
            has_cell_data = False

            if cds:
                has_point_data = True
                for cd in cds:
                    data = np.fromfile(cd, dtype=extension_to_dtype(extension(cd)))
                    name = os.path.basename(cd)
                    name = os.path.splitext(os.path.splitext(name)[0])[0]
                    name = name.replace(".", "_")
                    name_to_point_data[name] = data

                    if len(data) != len(data):
                        frame = inspect.currentframe()
                        print(
                            f"Error in {__file__} at line {frame.f_lineno}:\n .... data length is different from number of nodes {len(data)} != {len(data)}"
                        )
                        exit(1)

            cds = cell_data_steps[t]

            if cds:
                has_cell_data = True
                for cd in cds:
                    data = np.fromfile(cd, dtype=extension_to_dtype(extension(cd)))
                    name = os.path.basename(cd)
                    name = os.path.splitext(os.path.splitext(name)[0])[0]
                    name = name.replace(".", "_")
                    name_to_cell_data[name] = data

                    if len(data) != len(data):
                        frame = inspect.currentframe()
                        print(
                            f"Error in {__file__} at line {frame.f_lineno}:\n .... data length is different from number of nodes {len(data)} != {len(data)}"
                        )
                        exit(1)

            if has_point_data and not has_cell_data:
                writer.write_data(time_whole[t], point_data=name_to_point_data)
            elif not has_point_data and has_cell_data:
                writer.write_data(time_whole[t], cell_data=name_to_cell_data)
            elif has_point_data and has_cell_data:
                writer.write_data(
                    time_whole[t],
                    point_data=name_to_point_data,
                    cell_data=name_to_cell_data,
                )


def add_fields(field_data, storage, check_len):
    if field_data:
        paths = field_data.split(",")

        n_paths = len(paths)
        for i in range(0, n_paths):
            p = paths[i]
            files = glob.glob(p, recursive=False)
            files.sort()

            for f in files:
                field = Field(f)
                field.check_len(check_len)
                storage[field.name] = field.data


def raw_to_db(argv):
    usage = f"usage: {argv[0]} <input_folder> <output_mesh>"

    if len(argv) < 3:
        print(usage)
        sys.exit(1)

    raw_mesh_folder = argv[1]
    output_path = argv[2]
    raw_xyz_folder = raw_mesh_folder

    point_data = None
    cell_data = None

    transient = False
    time_step_format = "%s.%d.%d.raw"
    n_time_steps = 1
    time_whole = []

    cell_type = None
    verbose = False
    ssref = 0

    try:
        opts, args = getopt.getopt(
            argv[3:],
            "p:d:c:t:hv",
            [
                "coords=",
                "point_data=",
                "cell_type=",
                "cell_data=",
                "transient",
                "time_step_format",
                "n_time_steps=",
                "time_whole=",
                "time_whole_txt=",
                "help",
                "verbose",
                "ssref",
            ],
        )

    except getopt.GetoptError as err:
        print(err)
        print(usage)
        sys.exit(1)

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(usage)
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-p", "--point_data"):
            point_data = arg
        elif opt in ("-c", "--cell_data"):
            cell_data = arg
        elif opt in ("--transient"):
            transient = True
        elif opt in ("--time_whole"):
            time_whole = np.fromfile(arg, dtype=np.float32)
        elif opt in ("--time_whole_txt"):
            time_whole = np.loadtxt(arg, dtype=np.float32)
        elif opt in ("--time_step_format"):
            time_step_format = arg
        elif opt in ("--n_time_steps"):
            n_time_steps = int(arg)
        elif opt in ("--cell_type"):
            cell_type = arg
        elif opt in ("--coords"):
            raw_xyz_folder = arg
            if verbose:
                print(f"Using coords={arg}")
        elif opt in ("--ssref"):
            ssref = int(arg)

    mesh_meta_path = os.path.join(raw_mesh_folder, "meta.yaml")
    meta_type, mesh_meta, mesh_meta_path = element_type_from_meta(raw_mesh_folder)
    if cell_type is None:
        if meta_type is not None:
            cell_type = meta_type
            if verbose:
                print(f"Using element_type={cell_type} from {mesh_meta_path}")

    if transient:
        if len(time_whole) == 0:
            assert n_time_steps != 0
            time_whole = np.arange(0, n_time_steps)
        else:
            n_time_steps = len(time_whole)

        print(f"Found {n_time_steps} time steps!")

    points = []
    for pfn in ["x", "y", "z"]:
        path = detect_files(
            f"{raw_xyz_folder}/{pfn}.*", ["float16", "float32", "float64"]
        )

        if len(path) == 0:
            break

        path = path[0]
        if os.path.exists(path):
            if verbose:
                print(f"Reading {path}")
            x = np.fromfile(path, dtype=geom_t)
            points.append(x)

    # Attempt format x0, x1, x2
    if len(points) == 0:
        for d in range(0, 3):
            path = detect_files(
                f"{raw_xyz_folder}/x{d}.*", ["float16", "float32", "float64"]
            )
            if len(path) == 0:
                break

            path = path[0]
            dtype = extension_to_dtype(extension(path))
            if os.path.exists(path):
                if verbose:
                    print(f"Reading {path}")
                x = np.fromfile(path, dtype=dtype)
                points.append(x)

    idx = []
    for i in range(0, max_nodes_x_element):
        path = detect_files(
            f"{raw_mesh_folder}/i{i}.*", ["raw", "int16", "int32", "int64"]
        )
        if len(path) == 0:
            break

        path = path[0]
        dtype = extension_to_dtype(extension(path))
        if os.path.exists(path):
            if verbose:
                print(f"Reading {path}")
            ii = np.fromfile(path, dtype=dtype)
            idx.append(ii)

    if mesh_meta and "elem_num_nodes" in mesh_meta and len(idx) > 0:
        expected_nnodes = int(mesh_meta["elem_num_nodes"])
        if len(idx) != expected_nnodes:
            raise RuntimeError(
                f"meta.yaml elem_num_nodes={expected_nnodes} "
                f"but found {len(idx)} connectivity streams"
            )

    if cell_type in quad4_names:
        cell_type = "quad"

    if cell_type in quad9_names:
        cell_type = "quad9"

    if cell_type in hex8_names:
        cell_type = "hexahedron"

    reorder = None
    if cell_type in hex27_names:
        cell_type = "hexahedron27"
    elif cell_type in proteus_hex27_names:
        cell_type = "hexahedron27"
        reorder = proteus_hex27_to_hexahedron27

    if cell_type in tet4_names:
        cell_type = "tetra"

    if cell_type in tet10_names:
        cell_type = "tetra10"

    if cell_type in tri3_names:
        cell_type = "triangle"

    if cell_type in tri6_names:
        cell_type = "triangle6"

    # Do I need to do that?
    # if ssref > 1:
    #     # Convert ssmesh to standard mesh or to high-order rep
    #     assert cell_type != None

    #     if cell_type == "quad":
    #         idx, points = ssquad4_to_standard(ssref, idx, points)
    #     elif cell_type == "hexahedron"
    #         # Implement me!
    #         assert False

    if cell_type == None:
        if len(idx) == 3:
            cell_type = "triangle"
        elif len(idx) == 6:
            cell_type = "triangle6"
        elif len(idx) == 9:
            cell_type = "quad9"
        elif len(idx) == 4:
            if len(points) == 2:
                cell_type = "quad"
            else:
                cell_type = "tetra"
        elif len(idx) == 8:
            cell_type = "hexahedron"
        elif len(idx) == 27:
            cell_type = "hexahedron27"
        elif len(idx) == 10:
            cell_type = "tetra10"
        elif len(idx) == 2:
            cell_type = "line"
        elif len(idx) == 1:
            cell_type = "vertex"

    if cell_type == "quad" and len(idx) == 9:
        cell_type = "quad9"
    elif cell_type == "triangle" and len(idx) == 6:
        cell_type = "triangle6"
    elif cell_type == "tetra" and len(idx) == 10:
        cell_type = "tetra10"
    elif cell_type == "hexahedron" and len(idx) == 27:
        cell_type = "hexahedron27"

    print(f"numnodes = {len(idx)} -> {cell_type}")
    n_points = len(points[0])
    n_cells = len(idx[0])

    if n_points == 0 or n_cells == 0:
        print(f"Warning empty database at {raw_mesh_folder}")
        return

    points = np.array(points).transpose()
    cell_indices = np.array(idx).transpose()
    if reorder is not None:
        cell_indices = cell_indices[:, reorder]
    if cell_type == "hexahedron27":
        cell_indices = cell_indices[:, exodus_hex27_to_vtk_hex27]
        if verbose:
            print("Applied Exodus HEX27 -> VTK face-center node permutation")
    cells = [(cell_type, cell_indices)]

    if transient:
        print("Transient mode!")

        write_transient_data(
            output_path,
            points,
            cells,
            point_data,
            cell_data,
            n_time_steps,
            time_whole,
            time_step_format,
        )
    else:
        mesh = meshio.Mesh(points, cells)

        add_fields(point_data, mesh.point_data, n_points)
        add_fields(cell_data, mesh.cell_data, n_cells)

        if is_exodus_output(output_path) and verbose:
            print(
                "Writing Exodus via meshio (same VTK HEX27 layout as .vtu); "
                "use raw_to_exodusII for FEM/IOSS PATRAN ordering"
            )

        mesh.write(output_path)


# Example usage
# ./raw_to_db.py raw_db mesh_db.vtk --point_data='raw_db/point_data/*,raw_db/x.raw' --point_data_type='float64,float32'
if __name__ == "__main__":
    raw_to_db(sys.argv)
