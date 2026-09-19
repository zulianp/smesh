#!/usr/bin/env python3

import getopt
import os
import sys
from collections import OrderedDict

import numpy as np

from common.hex27_ordering import prepare_exodus_hex27_connectivity
from common.raw_io import (
    DEFAULT_GEOM_DTYPE,
    is_dtype_token,
    load_raw_mesh,
    read_array,
    sidesets_to_global,
    smesh_element_type,
    unique_preserve_order,
)
from common.raw_io import ELEMENT_INFO, exodus_element_type


try:
    geom_t
except NameError:
    print("raw_to_exodusII: self contained mode")
    geom_t = np.float32
    idx_t = np.int32
    element_idx_t = np.int32


LEN_STRING = 33
LEN_LINE = 81


def strip_typed_suffix(filename):
    parts = filename.split(".")
    if len(parts) >= 3 and parts[-1] == "raw" and is_dtype_token(parts[-2]):
        return ".".join(parts[:-2])
    if len(parts) >= 2 and is_dtype_token(parts[-1]):
        return ".".join(parts[:-1])
    if len(parts) >= 2 and parts[-1] == "raw":
        return ".".join(parts[:-1])
    return filename


def load_time_whole(folder):
    if not os.path.isdir(folder):
        return None
    for entry in sorted(os.listdir(folder)):
        if not entry.startswith("time_whole."):
            continue
        path = os.path.join(folder, entry)
        data = read_array(path, default_dtype=np.float32).astype(np.float32)
        if len(data) == 0:
            return None
        return data
    return None


def load_field_series(folder, expected_len):
    if not os.path.isdir(folder):
        return OrderedDict()

    fields = OrderedDict()
    for entry in sorted(os.listdir(folder)):
        path = os.path.join(folder, entry)
        if not os.path.isfile(path):
            continue

        stem = strip_typed_suffix(entry)
        parts = stem.split(".")
        time_index = 0
        if len(parts) >= 2 and parts[-1].isdigit():
            time_index = int(parts[-1])
            name = ".".join(parts[:-1])
        else:
            name = stem

        data = read_array(path)
        if len(data) != expected_len:
            raise RuntimeError(
                f"field {path} has invalid length {len(data)} (expected {expected_len})"
            )

        if name not in fields:
            fields[name] = {}
        fields[name][time_index] = data

    return fields


def align_time_series(point_fields, cell_fields, time_whole):
    max_time = 0
    for groups in (point_fields, cell_fields):
        for values in groups.values():
            if values:
                max_time = max(max_time, max(values.keys()))

    n_time_steps = max_time + 1 if max_time > 0 else 1

    if time_whole is None:
        time_whole = np.arange(n_time_steps, dtype=np.float32)
    else:
        time_whole = np.asarray(time_whole, dtype=np.float32)
        if len(time_whole) < n_time_steps:
            raise RuntimeError(
                f"time_whole has {len(time_whole)} entries but fields require {n_time_steps}"
            )
        n_time_steps = len(time_whole)

    for groups in (point_fields, cell_fields):
        for name, values in groups.items():
            for time_index in range(n_time_steps):
                if time_index not in values:
                    raise RuntimeError(f"missing time step {time_index} for field '{name}'")

    return time_whole, n_time_steps


def string_matrix(values, string_len=LEN_STRING):
    import netCDF4

    encoded = np.asarray([str(value)[:string_len] for value in values], dtype=f"S{string_len}")
    return netCDF4.stringtochar(encoded)


def _prop_ids(items, count):
    ids = []
    for i, item in enumerate(items):
        if isinstance(item, dict):
            exo_id = item.get("exodus_id")
        else:
            exo_id = getattr(item, "exodus_id", None)
        ids.append(int(exo_id) if exo_id is not None else i + 1)
    while len(ids) < count:
        ids.append(len(ids) + 1)
    return np.asarray(ids, dtype=np.int32)


def write_exodus(
    output_mesh,
    title,
    points,
    blocks,
    sidesets,
    nodesets,
    point_fields,
    cell_fields,
    time_whole,
):
    import netCDF4

    n_dim, n_nodes = points.shape
    n_elem = int(sum(block.connectivity.shape[0] for block in blocks))
    n_blocks = len(blocks)

    time_whole = np.asarray(time_whole, dtype=np.float32)
    n_time_steps = len(time_whole)

    point_field_names = list(point_fields.keys())
    cell_field_names = list(cell_fields.keys())
    coord_dtype = points.dtype if points.size else DEFAULT_GEOM_DTYPE

    begins = []
    offset = 0
    for block in blocks:
        begins.append(offset)
        offset += int(block.connectivity.shape[0])

    with netCDF4.Dataset(output_mesh, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.title = title
        nc.version = np.float32(5.1)
        nc.api_version = np.float32(5.1)
        nc.floating_point_word_size = np.int64(8)

        nc.createDimension("num_dim", n_dim)
        nc.createDimension("num_nodes", n_nodes)
        nc.createDimension("num_elem", n_elem)
        nc.createDimension("num_el_blk", max(n_blocks, 0))
        nc.createDimension("len_string", LEN_STRING)
        nc.createDimension("len_line", LEN_LINE)
        nc.createDimension("four", 4)
        nc.createDimension("time_step", None)

        if sidesets:
            nc.createDimension("num_side_sets", len(sidesets))
        if nodesets:
            nc.createDimension("num_node_sets", len(nodesets))

        time_var = nc.createVariable("time_whole", np.float32, ("time_step",))
        time_var[:] = time_whole

        coor_names = nc.createVariable("coor_names", "S1", ("num_dim", "len_string"))
        coor_names[:, :] = string_matrix(["x", "y", "z"][:n_dim])

        coord = nc.createVariable("coord", coord_dtype, ("num_dim", "num_nodes"))
        coord[:, :] = points

        if n_blocks:
            eb_status = nc.createVariable("eb_status", np.int32, ("num_el_blk",))
            eb_status[:] = np.ones(n_blocks, dtype=np.int32)

            eb_prop1 = nc.createVariable("eb_prop1", np.int32, ("num_el_blk",))
            eb_prop1.setncattr("name", "ID")
            eb_prop1[:] = _prop_ids(blocks, n_blocks)

            eb_names = nc.createVariable("eb_names", "S1", ("num_el_blk", "len_string"))
            eb_names[:, :] = string_matrix([block.name for block in blocks])

        for block_index, block in enumerate(blocks, start=1):
            conn = np.asarray(block.connectivity)
            smesh_type = smesh_element_type(block.element_type, nnodes=conn.shape[1])
            exo_type = exodus_element_type(smesh_type)
            info = ELEMENT_INFO[smesh_type]
            if conn.shape[1] != info["nnodes"]:
                raise RuntimeError(
                    f"connectivity width {conn.shape[1]} does not match {smesh_type} ({info['nnodes']})"
                )
            block_size = int(conn.shape[0])
            nnodes_per_elem = int(conn.shape[1])
            nc.createDimension(f"num_el_in_blk{block_index}", block_size)
            nc.createDimension(f"num_nod_per_el{block_index}", nnodes_per_elem)
            connect = nc.createVariable(
                f"connect{block_index}",
                np.int32,
                (f"num_el_in_blk{block_index}", f"num_nod_per_el{block_index}"),
            )
            connect.setncattr("elem_type", exo_type)
            if block_size:
                connect[:, :] = (conn.astype(np.int32, copy=False) + 1)

        if point_field_names:
            nc.createDimension("num_nod_var", len(point_field_names))
            name_nod_var = nc.createVariable("name_nod_var", "S1", ("num_nod_var", "len_string"))
            name_nod_var[:, :] = string_matrix(point_field_names)
            for index, name in enumerate(point_field_names, start=1):
                first = point_fields[name][0]
                values = nc.createVariable(
                    f"vals_nod_var{index}",
                    first.dtype,
                    ("time_step", "num_nodes"),
                )
                for time_index in range(n_time_steps):
                    values[time_index, :] = point_fields[name][time_index]

        if cell_field_names:
            nc.createDimension("num_elem_var", len(cell_field_names))
            name_elem_var = nc.createVariable(
                "name_elem_var",
                "S1",
                ("num_elem_var", "len_string"),
            )
            name_elem_var[:, :] = string_matrix(cell_field_names)

            elem_var_tab = nc.createVariable(
                "elem_var_tab",
                np.int32,
                ("num_el_blk", "num_elem_var"),
            )
            elem_var_tab[:, :] = np.ones((n_blocks, len(cell_field_names)), dtype=np.int32)

            for var_index, name in enumerate(cell_field_names, start=1):
                first = cell_fields[name][0]
                for block_index, block in enumerate(blocks, start=1):
                    begin = begins[block_index - 1]
                    end = begin + int(block.connectivity.shape[0])
                    values = nc.createVariable(
                        f"vals_elem_var{var_index}eb{block_index}",
                        first.dtype,
                        ("time_step", f"num_el_in_blk{block_index}"),
                    )
                    for time_index in range(n_time_steps):
                        values[time_index, :] = cell_fields[name][time_index][begin:end]

        if sidesets:
            ss_status = nc.createVariable("ss_status", np.int32, ("num_side_sets",))
            ss_status[:] = np.ones(len(sidesets), dtype=np.int32)

            ss_prop1 = nc.createVariable("ss_prop1", np.int32, ("num_side_sets",))
            ss_prop1.setncattr("name", "ID")
            ss_prop1[:] = _prop_ids(sidesets, len(sidesets))

            ss_names = nc.createVariable("ss_names", "S1", ("num_side_sets", "len_string"))
            ss_names[:, :] = string_matrix([sideset["name"] for sideset in sidesets])

            for ss_index, sideset in enumerate(sidesets, start=1):
                size = len(sideset["parent"])
                if len(sideset["lfi"]) != size:
                    raise RuntimeError(f"sideset '{sideset['name']}' has inconsistent lengths")
                ss_status[ss_index - 1] = 1 if size > 0 else 0
                # netCDF3 treats size 0 as NC_UNLIMITED (already used by time_step).
                if size == 0:
                    continue

                nc.createDimension(f"num_side_ss{ss_index}", size)
                elem_ss = nc.createVariable(
                    f"elem_ss{ss_index}",
                    np.int32,
                    (f"num_side_ss{ss_index}",),
                )
                side_ss = nc.createVariable(
                    f"side_ss{ss_index}",
                    np.int32,
                    (f"num_side_ss{ss_index}",),
                )
                elem_ss[:] = np.asarray(sideset["parent"], dtype=np.int64) + 1
                side_ss[:] = np.asarray(sideset["lfi"], dtype=np.int16) + 1

        if nodesets:
            ns_status = nc.createVariable("ns_status", np.int32, ("num_node_sets",))
            ns_status[:] = np.ones(len(nodesets), dtype=np.int32)

            ns_prop1 = nc.createVariable("ns_prop1", np.int32, ("num_node_sets",))
            ns_prop1.setncattr("name", "ID")
            ns_prop1[:] = _prop_ids(nodesets, len(nodesets))

            ns_names = nc.createVariable("ns_names", "S1", ("num_node_sets", "len_string"))
            ns_names[:, :] = string_matrix([nodeset.name for nodeset in nodesets])

            for ns_index, nodeset in enumerate(nodesets, start=1):
                nodes = unique_preserve_order(np.asarray(nodeset.nodes, dtype=np.int64))
                ns_status[ns_index - 1] = 1 if nodes.size > 0 else 0
                if nodes.size == 0:
                    continue
                nc.createDimension(f"num_nod_ns{ns_index}", len(nodes))
                node_ns = nc.createVariable(
                    f"node_ns{ns_index}",
                    np.int32,
                    (f"num_nod_ns{ns_index}",),
                )
                node_ns[:] = nodes + 1


def raw_to_exodusII(input_folder, output_mesh, title=None):
    mesh = load_raw_mesh(input_folder)
    for block in mesh.blocks:
        conn, element_type = prepare_exodus_hex27_connectivity(
            block.connectivity, block.element_type
        )
        block.connectivity = conn
        block.element_type = element_type

    offsets = mesh.element_offsets()
    sidesets = sidesets_to_global(mesh.sidesets, offsets)
    sidesets.sort(
        key=lambda item: (
            item["exodus_id"] is None,
            int(item["exodus_id"]) if item["exodus_id"] is not None else 0,
        )
    )
    nodesets = sorted(
        mesh.nodesets,
        key=lambda item: (
            item.exodus_id is None,
            int(item.exodus_id) if item.exodus_id is not None else 0,
        ),
    )

    point_fields = load_field_series(os.path.join(input_folder, "point_data"), mesh.n_nodes)
    cell_fields = load_field_series(os.path.join(input_folder, "cell_data"), mesh.n_elements)
    time_whole = load_time_whole(input_folder)
    time_whole, _ = align_time_series(point_fields, cell_fields, time_whole)

    if title is None:
        title = f"Created by raw_to_exodusII.py from {os.path.basename(os.path.abspath(input_folder))}"

    write_exodus(
        output_mesh=output_mesh,
        title=title,
        points=mesh.points,
        blocks=mesh.blocks,
        sidesets=sidesets,
        nodesets=nodesets,
        point_fields=point_fields,
        cell_fields=cell_fields,
        time_whole=time_whole,
    )


def main(argv):
    usage = f"usage: {argv[0]} <input_folder> <output_mesh> [--title=TITLE]"

    if len(argv) < 3:
        print(usage)
        sys.exit(1)

    input_folder = argv[1]
    output_mesh = argv[2]
    title = None

    try:
        opts, args = getopt.getopt(argv[3:], "h", ["title=", "help"])
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
            sys.exit(0)
        if opt == "--title":
            title = arg

    raw_to_exodusII(input_folder, output_mesh, title=title)


if __name__ == "__main__":
    main(sys.argv)
