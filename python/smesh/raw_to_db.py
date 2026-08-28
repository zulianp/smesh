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
wedge6_names = ("wedge", "wedge6", "WEDGE6", "WEDGE", "prism", "prism6")
pyramid5_names = ("pyramid", "pyramid5", "PYRAMID5", "PYRAMID")

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


def parse_smesh_block_list(meta_path):
    """Parse per-block name/type from smesh multi-block meta.yaml."""
    blocks = []
    if not os.path.exists(meta_path):
        return blocks

    current = None
    in_blocks = False
    with open(meta_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped == "blocks:":
                in_blocks = True
                continue
            if not in_blocks:
                continue
            if stripped.startswith("- name:"):
                if current:
                    blocks.append(current)
                current = {"name": stripped.split(":", 1)[1].strip()}
                continue
            if current is None:
                continue
            if (
                not line.startswith(" ")
                and not line.startswith("\t")
                and not stripped.startswith("-")
                and ":" in stripped
            ):
                blocks.append(current)
                current = None
                in_blocks = False
                continue
            if ":" not in stripped or stripped.startswith("- "):
                continue
            key, value = stripped.split(":", 1)
            key = key.strip()
            value = value.strip()
            if key in ("name", "element_type", "cell_type", "elem_type"):
                current[key] = value
            elif key in ("elem_num_nodes", "n_elements"):
                current[key] = int(value)
    if current:
        blocks.append(current)
    return blocks


def discover_block_dirs(mesh_folder):
    blocks_root = os.path.join(mesh_folder, "blocks")
    if not os.path.isdir(blocks_root):
        return []
    names = []
    for name in sorted(os.listdir(blocks_root)):
        folder = os.path.join(blocks_root, name)
        if not os.path.isdir(folder):
            continue
        if detect_files(f"{folder}/i0.*", ["raw", "int16", "int32", "int64"]):
            names.append(name)
    return names


def load_idx_from_folder(folder, verbose=False):
    idx = []
    for i in range(0, max_nodes_x_element):
        path = detect_files(
            f"{folder}/i{i}.*", ["raw", "int16", "int32", "int64"]
        )
        if len(path) == 0:
            break
        path = path[0]
        dtype = extension_to_dtype(extension(path))
        if os.path.exists(path):
            if verbose:
                print(f"Reading {path}")
            idx.append(np.fromfile(path, dtype=dtype))
    return idx


def resolve_meshio_cell_type(cell_type, nnodes, n_spatial_dim):
    reorder = None
    if cell_type in quad4_names:
        cell_type = "quad"
    if cell_type in quad9_names:
        cell_type = "quad9"
    if cell_type in hex8_names:
        cell_type = "hexahedron"
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
    if cell_type in wedge6_names:
        cell_type = "wedge"
    if cell_type in pyramid5_names:
        cell_type = "pyramid"

    if cell_type is None:
        if nnodes == 3:
            cell_type = "triangle"
        elif nnodes == 6:
            cell_type = "triangle6"
        elif nnodes == 9:
            cell_type = "quad9"
        elif nnodes == 4:
            cell_type = "quad" if n_spatial_dim == 2 else "tetra"
        elif nnodes == 5:
            cell_type = "pyramid"
        elif nnodes == 8:
            cell_type = "hexahedron"
        elif nnodes == 27:
            cell_type = "hexahedron27"
        elif nnodes == 10:
            cell_type = "tetra10"
        elif nnodes == 2:
            cell_type = "line"
        elif nnodes == 1:
            cell_type = "vertex"

    if cell_type == "quad" and nnodes == 9:
        cell_type = "quad9"
    elif cell_type == "triangle" and nnodes == 6:
        cell_type = "triangle6"
    elif cell_type == "tetra" and nnodes == 10:
        cell_type = "tetra10"
    elif cell_type == "hexahedron" and nnodes == 27:
        cell_type = "hexahedron27"

    return cell_type, reorder


def connectivity_to_cell_block(idx, cell_type, n_spatial_dim, verbose=False):
    mio, reorder = resolve_meshio_cell_type(cell_type, len(idx), n_spatial_dim)
    if mio is None:
        raise RuntimeError(
            f"unable to infer cell type from {len(idx)} connectivity streams"
        )
    cell_indices = np.array(idx).transpose()
    if reorder is not None:
        cell_indices = cell_indices[:, reorder]
    if mio == "hexahedron27":
        cell_indices = cell_indices[:, exodus_hex27_to_vtk_hex27]
        if verbose:
            print("Applied Exodus HEX27 -> VTK face-center node permutation")
    return mio, cell_indices


def element_type_from_meta(folder):
    meta_path = os.path.join(folder, "meta.yaml")
    meta = read_simple_meta(meta_path)
    n_blocks = int(meta["n_blocks"]) if meta.get("n_blocks") else 0
    if n_blocks > 1:
        return None, meta, meta_path
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
    def read_transient_field(path, center_size):
        data = np.fromfile(path, dtype=extension_to_dtype(extension(path)))
        if ".vec3." in os.path.basename(path):
            if len(data) != 3 * center_size:
                print(
                    f"Invalid vec3 field length {len(data)} for {path}; expected {3 * center_size}"
                )
                sys.exit(1)
            data = data.reshape((center_size, 3))
        return data

    def transient_field_name(path):
        name = os.path.basename(path)
        name = os.path.splitext(os.path.splitext(name)[0])[0]
        name = name.replace(".vec3", "")
        return name.replace(".", "_")

    with meshio.xdmf.TimeSeriesWriter(output_path) as writer:
        writer.write_points_cells(points, cells)
        n_points = len(points)
        n_cells = sum(block[1].shape[0] for block in cells)
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
                    data = read_transient_field(cd, n_points)
                    name = transient_field_name(cd)
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
                    data = read_transient_field(cd, n_cells)
                    name = transient_field_name(cd)
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


def merge_cells_by_type(cell_blocks, cell_data):
    """One connectivity array per VTK cell type so ParaView attaches CELL_DATA."""
    types = []
    conn_by_type = {}
    index_by_type = {}
    for i, block in enumerate(cell_blocks):
        cell_type = block.type if hasattr(block, "type") else block[0]
        conn = block.data if hasattr(block, "data") else block[1]
        if cell_type not in conn_by_type:
            types.append(cell_type)
            conn_by_type[cell_type] = []
            index_by_type[cell_type] = []
        conn_by_type[cell_type].append(np.asarray(conn))
        index_by_type[cell_type].append(i)

    new_cells = [(t, np.concatenate(conn_by_type[t], axis=0)) for t in types]
    new_cell_data = {}
    for name, pieces in cell_data.items():
        seq = pieces if isinstance(pieces, (list, tuple)) else [pieces]
        merged = []
        for t in types:
            parts = [np.asarray(seq[i]) for i in index_by_type[t] if i < len(seq)]
            merged.append(np.concatenate(parts, axis=0))
        new_cell_data[name] = merged
    return new_cells, new_cell_data


def write_paraview_mesh(mesh, output_path):
    # VTK 5.1 FIELD + OFFSETS often hides cell arrays in ParaView.
    ext = os.path.splitext(output_path)[1].lower()
    if ext == ".vtk":
        import meshio.vtk as meshio_vtk

        meshio_vtk.write(output_path, mesh, fmt_version="4.2")
    else:
        mesh.write(output_path)


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

    n_spatial_dim = len(points)
    block_specs = parse_smesh_block_list(mesh_meta_path)
    spec_by_name = {spec["name"]: spec for spec in block_specs if "name" in spec}
    block_names = discover_block_dirs(raw_mesh_folder)

    cells = []
    block_sizes = []
    if block_names:
        cli_cell_type = cell_type
        for name in block_names:
            folder = os.path.join(raw_mesh_folder, "blocks", name)
            idx = load_idx_from_folder(folder, verbose)
            if not idx:
                print(f"Error: no connectivity in {folder}")
                sys.exit(1)
            spec = spec_by_name.get(name, {})
            et = spec.get("element_type") or spec.get("cell_type") or spec.get(
                "elem_type"
            )
            if len(block_names) == 1 and cli_cell_type is not None:
                et = cli_cell_type
            if "elem_num_nodes" in spec and int(spec["elem_num_nodes"]) != len(idx):
                raise RuntimeError(
                    f"block '{name}' elem_num_nodes={spec['elem_num_nodes']} "
                    f"but found {len(idx)} connectivity streams"
                )
            mio, cell_indices = connectivity_to_cell_block(
                idx, et, n_spatial_dim, verbose
            )
            print(f"block '{name}': numnodes = {len(idx)} -> {mio}")
            cells.append((mio, cell_indices))
            block_sizes.append(cell_indices.shape[0])
    else:
        idx = load_idx_from_folder(raw_mesh_folder, verbose)
        if mesh_meta and "elem_num_nodes" in mesh_meta and len(idx) > 0:
            expected_nnodes = int(mesh_meta["elem_num_nodes"])
            if len(idx) != expected_nnodes:
                raise RuntimeError(
                    f"meta.yaml elem_num_nodes={expected_nnodes} "
                    f"but found {len(idx)} connectivity streams"
                )
        if not idx:
            print(f"Error: no connectivity in {raw_mesh_folder}")
            sys.exit(1)
        mio, cell_indices = connectivity_to_cell_block(
            idx, cell_type, n_spatial_dim, verbose
        )
        print(f"numnodes = {len(idx)} -> {mio}")
        cells.append((mio, cell_indices))
        block_sizes.append(cell_indices.shape[0])

    n_points = len(points[0])
    n_cells = int(sum(block_sizes))

    if n_points == 0 or n_cells == 0:
        print(f"Warning empty database at {raw_mesh_folder}")
        return

    points = np.array(points).transpose()

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
        if len(block_sizes) == 1:
            add_fields(cell_data, mesh.cell_data, n_cells)
        elif cell_data:
            raw_cell = {}
            add_fields(cell_data, raw_cell, n_cells)
            for name, data in raw_cell.items():
                split = []
                offset = 0
                for sz in block_sizes:
                    split.append(data[offset : offset + sz])
                    offset += sz
                mesh.cell_data[name] = split
        if len(block_sizes) > 1:
            mesh.cell_data["block"] = [
                np.full(sz, i, dtype=np.int32) for i, sz in enumerate(block_sizes)
            ]

        if is_exodus_output(output_path) and verbose:
            print(
                "Writing Exodus via meshio (same VTK HEX27 layout as .vtu); "
                "use raw_to_exodusII for FEM/IOSS PATRAN ordering"
            )

        cells, cell_data_out = merge_cells_by_type(mesh.cells, mesh.cell_data)
        mesh = meshio.Mesh(
            mesh.points,
            cells,
            point_data=mesh.point_data,
            cell_data=cell_data_out,
        )
        write_paraview_mesh(mesh, output_path)


# Example usage
# ./raw_to_db.py raw_db mesh_db.vtk --point_data='raw_db/point_data/*,raw_db/x.raw' --point_data_type='float64,float32'
if __name__ == "__main__":
    raw_to_db(sys.argv)
