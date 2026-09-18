#!/usr/bin/env python3

import meshio
import numpy as np
import sys, getopt
import os

from common.hex27_ordering import vtk_hex27_to_exodus_hex27
from common.raw_io import (
    DEFAULT_IDX_DTYPE,
    Block,
    Nodeset,
    RawMesh,
    is_exodus_path,
    mkdir,
    smesh_type_from_meshio,
    write_raw_mesh,
)
from common.utils import dtype_to_extension


try:
    geom_t
except NameError:
    geom_t = None
    idx_t = np.int32
    real_t = np.float64


def flatten_cell_data(data, cells, elem_type_filter):
    """meshio cell_data is a list of arrays, one per cell block."""
    if data is None:
        return None
    if isinstance(data, np.ndarray) and data.dtype != object:
        return np.asarray(data)

    seq = data if isinstance(data, (list, tuple)) else [data]
    parts = []
    for i, block in enumerate(cells):
        if elem_type_filter is not None and block.type != elem_type_filter:
            continue
        if i >= len(seq):
            break
        parts.append(np.asarray(seq[i]))
    if not parts:
        return None
    return np.concatenate(parts, axis=0)


def _cell_block_conn(block):
    return np.asarray(block.data if hasattr(block, "data") else block[1])


def _cell_block_type(block):
    return block.type if hasattr(block, "type") else block[0]


def block_names_from_meshio(mesh, selected):
    n_blocks = len(selected)
    names = [None] * n_blocks
    cell_sets = getattr(mesh, "cell_sets", None) or {}
    for name, pieces in cell_sets.items():
        for local_i, (global_i, block) in enumerate(selected):
            if global_i >= len(pieces):
                continue
            idx = np.asarray(pieces[global_i])
            ncells = _cell_block_conn(block).shape[0]
            if idx.size == ncells and (
                ncells == 0 or (idx.min() == 0 and idx.max() == ncells - 1)
            ):
                if names[local_i] is None:
                    names[local_i] = str(name)
    used = set()
    out = []
    for i, name in enumerate(names):
        if not name:
            name = f"block_{i}"
        original = name
        suffix = 1
        while name in used:
            suffix += 1
            name = f"{original}_{suffix}"
        used.add(name)
        out.append(name)
    return out


def db_to_raw(argv):
    usage = f"usage: {argv[0]} <input_mesh> <output_folder>"
    if len(argv) < 3:
        print(usage)
        sys.exit(1)

    input_mesh = argv[1]
    output_folder = argv[2]
    elem_type_filter = None

    try:
        opts, args = getopt.getopt(argv[3:], "e:h", ["select_elem_type=", "help"])
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
        elif opt in ("-e", "--select_elem_type"):
            elem_type_filter = arg

    if is_exodus_path(input_mesh):
        from exodusII_to_raw import exodusII_to_raw

        exodusII_to_raw(input_mesh, output_folder)
        return

    mkdir(output_folder)
    mesh = meshio.read(input_mesh)

    selected = []
    for i, block in enumerate(mesh.cells):
        if elem_type_filter is not None and _cell_block_type(block) != elem_type_filter:
            continue
        selected.append((i, block))

    if not selected:
        print(f"no cell blocks selected from {input_mesh}")
        sys.exit(1)

    spatial_dim = int(mesh.points.shape[1])
    names = block_names_from_meshio(mesh, selected)
    blocks = []
    idx_dtype = DEFAULT_IDX_DTYPE
    for name, (_i, cell_block) in zip(names, selected):
        conn = _cell_block_conn(cell_block).copy()
        cell_type = _cell_block_type(cell_block)
        if cell_type == "hexahedron27":
            conn = conn[:, vtk_hex27_to_exodus_hex27]
        smesh_type = smesh_type_from_meshio(
            cell_type, nnodes=conn.shape[1], spatial_dim=spatial_dim
        )
        blocks.append(
            Block(
                name=name,
                element_type=smesh_type,
                connectivity=conn.astype(idx_dtype, copy=False),
            )
        )
        print(f"block '{name}': {smesh_type} n_elements={conn.shape[0]}")

    geom_dtype = np.dtype(mesh.points.dtype if geom_t is None else geom_t)
    points = np.transpose(np.asarray(mesh.points, dtype=geom_dtype))

    nodesets = []
    point_sets = getattr(mesh, "point_sets", None) or {}
    for name, nodes in point_sets.items():
        nodesets.append(Nodeset(name=str(name), nodes=np.asarray(nodes, dtype=idx_dtype)))

    raw = RawMesh(
        points=points,
        blocks=blocks,
        nodesets=nodesets,
        idx_dtype=idx_dtype,
        geom_dtype=geom_dtype,
    )
    write_raw_mesh(output_folder, raw)

    point_data_dir = os.path.join(output_folder, "point_data")
    if mesh.point_data:
        mkdir(point_data_dir)
        for key in mesh.point_data:
            data = np.asarray(mesh.point_data[key])
            data.tofile(os.path.join(point_data_dir, f"{key}.{dtype_to_extension(data.dtype)}"))

    if mesh.cell_data:
        cell_data_dir = os.path.join(output_folder, "cell_data")
        mkdir(cell_data_dir)
        print("Cell data:")
        for key in mesh.cell_data:
            arr = flatten_cell_data(mesh.cell_data[key], mesh.cells, elem_type_filter)
            if arr is None:
                print(f"\t- {key} (unable to convert)")
                continue
            ext = dtype_to_extension(arr.dtype)
            path = os.path.join(cell_data_dir, f"{key}.{ext}")
            np.ascontiguousarray(arr).tofile(path)
            print(f"\t- {key} {arr.shape} {arr.dtype} -> {path}")


if __name__ == "__main__":
    db_to_raw(sys.argv)
