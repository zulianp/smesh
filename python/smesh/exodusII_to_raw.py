#!/usr/bin/env python3

import getopt
import os
import sys

import numpy as np

from common.raw_io import (
    DEFAULT_IDX_DTYPE,
    DEFAULT_PARENT_DTYPE,
    Block,
    Nodeset,
    RawMesh,
    element_offsets_from_counts,
    mkdir,
    smesh_element_type,
    split_global_sideset,
    write_raw_mesh,
)
from common.utils import dtype_to_extension


try:
    geom_t
except NameError:
    print("exodusII_to_raw: self contained mode")
    geom_t = None
    idx_t = np.int32
    element_idx_t = np.int32


def decode_name(value, fallback):
    if value is None:
        return fallback
    if isinstance(value, bytes):
        value = value.decode("ascii", "ignore")
    text = str(value).replace("\x00", "").strip()
    return text if text else fallback


def read_exodus_names(nc, var_name, count, fallback_prefix):
    names = []
    if var_name in nc.variables and count > 0:
        raw = nc.variables[var_name][:]
        try:
            import netCDF4

            decoded = netCDF4.chartostring(raw)
        except Exception:
            decoded = raw
        if count == 1 and np.ndim(decoded) == 0:
            decoded = [decoded]
        for i in range(count):
            item = decoded[i] if i < len(decoded) else None
            names.append(decode_name(item, f"{fallback_prefix}{i + 1}"))
    else:
        names = [f"{fallback_prefix}{i + 1}" for i in range(count)]
    return names


def read_prop_ids(nc, var_name, count):
    if var_name not in nc.variables or count <= 0:
        return [i + 1 for i in range(count)]
    values = np.asarray(nc.variables[var_name][:]).reshape(-1)
    ids = []
    for i in range(count):
        ids.append(int(values[i]) if i < values.size else i + 1)
    return ids


def load_exodus_coords(nc):
    if "coord" in nc.variables:
        coords = np.asarray(nc.variables["coord"][:])
        return coords
    axes = []
    for name in ("coordx", "coordy", "coordz"):
        if name in nc.variables:
            axes.append(np.asarray(nc.variables[name][:]))
    if not axes:
        raise RuntimeError("no coordinates in Exodus file")
    return np.vstack(axes)


def write_exodus_fields(nc, output_folder):
    n_time_steps = 1
    if "time_whole" in nc.variables:
        time_whole = np.asarray(nc.variables["time_whole"][:]).astype(np.float32, copy=False)
        n_time_steps = time_whole.shape[0] if time_whole.ndim else 1
        if time_whole.size:
            time_whole.tofile(
                os.path.join(output_folder, f"time_whole.{dtype_to_extension(np.float32)}")
            )
    print(f"n_time_steps = {n_time_steps}")

    if "name_nod_var" not in nc.variables:
        return
    import netCDF4

    name_nod_var = nc.variables["name_nod_var"]
    nvars = name_nod_var.shape[0]
    print(f"Point data, nvars = {nvars}")
    point_data_dir = os.path.join(output_folder, "point_data")
    mkdir(point_data_dir)
    for i in range(nvars):
        var_key = f"vals_nod_var{i + 1}"
        if var_key not in nc.variables:
            continue
        var = nc.variables[var_key]
        var_name = decode_name(netCDF4.chartostring(name_nod_var[i, :]), f"nod_var{i + 1}")
        print(f" - {var_name}, dtype {var.dtype}")
        prefix = os.path.join(point_data_dir, str(var_name))
        if n_time_steps <= 1:
            np.asarray(var[:]).tofile(f"{prefix}.{dtype_to_extension(var.dtype)}")
        else:
            pad = int(np.ceil(np.log10(max(n_time_steps, 1))))
            fmt = f"%s.%0.{pad}d.{dtype_to_extension(var.dtype)}"
            for t in range(n_time_steps):
                np.asarray(var[t, :]).tofile(fmt % (prefix, t))


def load_exodus_blocks(nc):
    num_el_blk = nc.dimensions["num_el_blk"].size if "num_el_blk" in nc.dimensions else 0
    names = read_exodus_names(nc, "eb_names", num_el_blk, "block_")
    ids = read_prop_ids(nc, "eb_prop1", num_el_blk)
    used_names = set()
    blocks = []
    idx_dtype = DEFAULT_IDX_DTYPE
    for b in range(num_el_blk):
        var_name = f"connect{b + 1}"
        if var_name not in nc.variables:
            raise RuntimeError(f"missing {var_name}")
        connect = nc.variables[var_name]
        conn = np.asarray(connect[:])
        if conn.ndim != 2:
            conn = conn.reshape((-1, conn.size if conn.size else 0))
        nnodes = int(conn.shape[1]) if conn.ndim == 2 else 0
        raw_type = getattr(connect, "elem_type", None)
        element_type = smesh_element_type(raw_type, nnodes=nnodes if nnodes else None)
        name = names[b] if b < len(names) else f"block_{b + 1}"
        if name in used_names:
            name = f"{name}_{b + 1}"
        used_names.add(name)
        zero_based = conn.astype(idx_dtype, copy=False) - np.array(1, dtype=idx_dtype)
        blocks.append(
            Block(
                name=name,
                element_type=element_type,
                connectivity=zero_based,
                exodus_id=ids[b] if b < len(ids) else b + 1,
            )
        )
        print(f"block '{name}': elem_type={element_type} n_elements={zero_based.shape[0]}")
    return blocks


def load_exodus_sidesets(nc, offsets):
    nss = nc.dimensions["num_side_sets"].size if "num_side_sets" in nc.dimensions else 0
    names = read_exodus_names(nc, "ss_names", nss, "sideset")
    ids = read_prop_ids(nc, "ss_prop1", nss)
    sidesets = []
    print(f"num_sidesets={nss}")
    for i in range(nss):
        name = names[i]
        ssidx = i + 1
        elem_key = f"elem_ss{ssidx}"
        side_key = f"side_ss{ssidx}"
        if elem_key not in nc.variables or side_key not in nc.variables:
            print(f"sideset = {name} size=0")
            sidesets.extend(
                split_global_sideset(
                    name,
                    np.zeros((0,), dtype=np.int64),
                    np.zeros((0,), dtype=np.int16),
                    offsets,
                    exodus_id=ids[i] if i < len(ids) else ssidx,
                )
            )
            continue
        elem = np.asarray(nc.variables[elem_key][:], dtype=np.int64) - 1
        side = np.asarray(nc.variables[side_key][:], dtype=np.int16) - 1
        print(f"sideset = {name} size={elem.size}")
        sidesets.extend(
            split_global_sideset(name, elem, side, offsets, exodus_id=ids[i] if i < len(ids) else ssidx)
        )
    return sidesets


def load_exodus_nodesets(nc):
    nns = nc.dimensions["num_node_sets"].size if "num_node_sets" in nc.dimensions else 0
    names = read_exodus_names(nc, "ns_names", nns, "nodeset")
    ids = read_prop_ids(nc, "ns_prop1", nns)
    nodesets = []
    print(f"num_nodesets={nns}")
    for i in range(nns):
        name = names[i]
        key = f"node_ns{i + 1}"
        if key not in nc.variables:
            nodes = np.zeros((0,), dtype=DEFAULT_IDX_DTYPE)
        else:
            nodes = np.asarray(nc.variables[key][:], dtype=DEFAULT_IDX_DTYPE) - 1
        print(f"nodeset = {name} size={nodes.size}")
        nodesets.append(
            Nodeset(
                name=name,
                nodes=nodes,
                exodus_id=ids[i] if i < len(ids) else i + 1,
            )
        )
    return nodesets


def exodusII_to_raw(input_mesh, output_folder):
    import netCDF4

    mkdir(output_folder)
    nc = netCDF4.Dataset(input_mesh)

    coords = load_exodus_coords(nc)
    geom_dtype = np.dtype(coords.dtype) if geom_t is None else np.dtype(geom_t)
    points = np.asarray(coords, dtype=geom_dtype)

    print(f"num_elem = {nc.dimensions['num_elem'].size if 'num_elem' in nc.dimensions else 0}")
    print(f"num_el_blk = {nc.dimensions['num_el_blk'].size if 'num_el_blk' in nc.dimensions else 0}")

    write_exodus_fields(nc, output_folder)

    blocks = load_exodus_blocks(nc)
    offsets = element_offsets_from_counts([block.n_elements for block in blocks])
    sidesets = load_exodus_sidesets(nc, offsets)
    nodesets = load_exodus_nodesets(nc)
    nc.close()

    mesh = RawMesh(
        points=points,
        blocks=blocks,
        sidesets=sidesets,
        nodesets=nodesets,
        idx_dtype=DEFAULT_IDX_DTYPE,
        geom_dtype=geom_dtype,
        parent_dtype=DEFAULT_PARENT_DTYPE,
    )
    write_raw_mesh(output_folder, mesh)


def main(argv):
    usage = f"usage: {argv[0]} <input_mesh> <output_folder>"
    if len(argv) < 3:
        print(usage)
        sys.exit(1)

    input_mesh = argv[1]
    output_folder = argv[2]
    try:
        opts, args = getopt.getopt(argv[3:], "h", ["help"])
    except getopt.GetoptError as err:
        print(err)
        print(usage)
        sys.exit(1)

    if args:
        print(f"unexpected positional arguments: {' '.join(args)}")
        print(usage)
        sys.exit(1)

    for opt, _arg in opts:
        if opt in ("-h", "--help"):
            print(usage)
            sys.exit(0)

    exodusII_to_raw(input_mesh, output_folder)


if __name__ == "__main__":
    main(sys.argv)
