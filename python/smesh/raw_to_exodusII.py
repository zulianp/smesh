#!/usr/bin/env python3

import getopt
import os
import re
import sys
from collections import OrderedDict

import numpy as np


try:
    geom_t
except NameError:
    print("raw_to_exodusII: self contained mode")
    geom_t = np.float32
    idx_t = np.int32
    element_idx_t = np.int32


LEN_STRING = 33
LEN_LINE = 81

ELEMENT_INFO = {
    "HEX8": {"exodus": "HEX8", "nnodes": 8},
    "HEX": {"exodus": "HEX8", "nnodes": 8},
    "TET4": {"exodus": "TETRA", "nnodes": 4},
    "TETRA": {"exodus": "TETRA", "nnodes": 4},
    "tetra": {"exodus": "TETRA", "nnodes": 4},
    "tetra4": {"exodus": "TETRA", "nnodes": 4},
    "QUAD4": {"exodus": "QUAD4", "nnodes": 4},
    "QUAD": {"exodus": "QUAD4", "nnodes": 4},
    "quad4": {"exodus": "QUAD4", "nnodes": 4},
    "TRI3": {"exodus": "TRI3", "nnodes": 3},
    "TRI": {"exodus": "TRI3", "nnodes": 3},
    "tri3": {"exodus": "TRI3", "nnodes": 3},
}


def is_dtype_token(token):
    try:
        np.dtype(token)
        return True
    except TypeError:
        return False


def strip_typed_suffix(filename):
    parts = filename.split(".")
    if len(parts) >= 3 and parts[-1] == "raw" and is_dtype_token(parts[-2]):
        return ".".join(parts[:-2])
    if len(parts) >= 2 and is_dtype_token(parts[-1]):
        return ".".join(parts[:-1])
    if len(parts) >= 2 and parts[-1] == "raw":
        return ".".join(parts[:-1])
    return filename


def dtype_from_path(path, default_dtype=None):
    name = os.path.basename(path)
    parts = name.split(".")
    if len(parts) >= 3 and parts[-1] == "raw" and is_dtype_token(parts[-2]):
        return np.dtype(parts[-2])
    if len(parts) >= 2 and is_dtype_token(parts[-1]):
        return np.dtype(parts[-1])
    return np.dtype(default_dtype) if default_dtype is not None else None


def read_array(path, default_dtype=None, count=None):
    dtype = dtype_from_path(path, default_dtype)
    if dtype is None:
        if count is not None:
            nbytes = os.path.getsize(path)
            if count == 0:
                dtype = np.dtype(default_dtype if default_dtype is not None else np.int32)
            else:
                itemsize = nbytes // count
                if itemsize == 2:
                    dtype = np.int16
                elif itemsize == 4:
                    dtype = np.int32
                elif itemsize == 8:
                    dtype = np.int64
                else:
                    raise RuntimeError(f"unable to infer dtype for {path}")
        else:
            dtype = np.dtype(default_dtype if default_dtype is not None else np.float32)
    return np.fromfile(path, dtype=dtype)


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


def normalize_element_type(raw_type):
    if raw_type is None:
        raise RuntimeError("missing element_type in meta.yaml")

    if raw_type not in ELEMENT_INFO:
        raise RuntimeError(f"unsupported element_type '{raw_type}'")

    return raw_type, ELEMENT_INFO[raw_type]


def find_axis_path(folder, axis):
    candidates = []
    for entry in os.listdir(folder):
        if not entry.startswith(f"{axis}."):
            continue
        if entry.endswith(".yaml"):
            continue
        candidates.append(os.path.join(folder, entry))

    candidates.sort()
    return candidates[0] if candidates else None


def load_points(folder):
    points = []
    for axis in ("x", "y", "z"):
        path = find_axis_path(folder, axis)
        if path is None:
            break
        points.append(read_array(path, default_dtype=geom_t).astype(geom_t))

    if not points:
        raise RuntimeError(f"no coordinate files found in {folder}")

    n_nodes = len(points[0])
    for axis in points:
        if len(axis) != n_nodes:
            raise RuntimeError("coordinate arrays have inconsistent lengths")

    return np.vstack(points)


def load_connectivity(folder):
    index_pattern = re.compile(r"^i(\d+)\.")
    entries = []
    for entry in os.listdir(folder):
        match = index_pattern.match(entry)
        if match:
            entries.append((int(match.group(1)), os.path.join(folder, entry)))

    if not entries:
        raise RuntimeError(f"no connectivity files found in {folder}")

    entries.sort(key=lambda item: item[0])
    arrays = [read_array(path, default_dtype=idx_t).astype(idx_t) for _, path in entries]

    n_elements = len(arrays[0])
    for array in arrays:
        if len(array) != n_elements:
            raise RuntimeError("connectivity arrays have inconsistent lengths")

    return np.column_stack(arrays)


def load_block_ranges(folder, n_elements):
    blocks_dir = os.path.join(folder, "blocks")
    if not os.path.isdir(blocks_dir):
        return [{"name": "block_1", "begin": 0, "end": n_elements}]

    blocks = []
    for entry in sorted(os.listdir(blocks_dir)):
        path = os.path.join(blocks_dir, entry)
        if not os.path.isfile(path):
            continue
        values = read_array(path, default_dtype=np.int64)
        if len(values) != 2:
            raise RuntimeError(f"invalid block range in {path}")
        begin = int(values[0])
        end = int(values[1])
        if begin < 0 or end < begin or end > n_elements:
            raise RuntimeError(f"invalid block bounds [{begin}, {end}) in {path}")
        blocks.append({"name": strip_typed_suffix(entry), "begin": begin, "end": end})

    if not blocks:
        return [{"name": "block_1", "begin": 0, "end": n_elements}]

    return blocks


def load_sidesets(folder):
    sidesets_dir = os.path.join(folder, "sidesets")
    if not os.path.isdir(sidesets_dir):
        return []

    sidesets = []
    for name in sorted(os.listdir(sidesets_dir)):
        ss_dir = os.path.join(sidesets_dir, name)
        if not os.path.isdir(ss_dir):
            continue

        lfi_path = None
        parent_path = None
        node_paths = []

        for entry in sorted(os.listdir(ss_dir)):
            path = os.path.join(ss_dir, entry)
            if not os.path.isfile(path):
                continue
            stem = strip_typed_suffix(entry)
            if stem == "lfi":
                lfi_path = path
            elif stem == "parent":
                parent_path = path
            else:
                prefix = f"{name}."
                if entry.startswith(prefix):
                    suffix = entry[len(prefix) :]
                    suffix = strip_typed_suffix(suffix)
                    if suffix.isdigit():
                        node_paths.append((int(suffix), path))

        if lfi_path is None or parent_path is None:
            raise RuntimeError(f"incomplete sideset data in {ss_dir}")

        lfi = read_array(lfi_path, default_dtype=np.int16).astype(np.int16)
        count = len(lfi)
        parent = read_array(parent_path, default_dtype=element_idx_t, count=count).astype(
            element_idx_t
        )

        node_paths.sort(key=lambda item: item[0])
        node_columns = []
        for _, path in node_paths:
            node_columns.append(read_array(path, default_dtype=idx_t, count=count).astype(idx_t))

        sidesets.append(
            {
                "name": name,
                "parent": parent,
                "lfi": lfi,
                "nodes": node_columns,
            }
        )

    return sidesets


def load_nodesets(folder, sidesets):
    nodesets_dir = os.path.join(folder, "nodesets")
    nodesets = []

    if os.path.isdir(nodesets_dir):
        for entry in sorted(os.listdir(nodesets_dir)):
            path = os.path.join(nodesets_dir, entry)
            if os.path.isdir(path):
                files = [
                    os.path.join(path, child)
                    for child in sorted(os.listdir(path))
                    if os.path.isfile(os.path.join(path, child))
                ]
                candidates = []
                for candidate in files:
                    stem = strip_typed_suffix(os.path.basename(candidate))
                    if stem in ("nodeset", "nodes", entry):
                        candidates.append(candidate)
                if not candidates:
                    candidates = files
                if not candidates:
                    continue
                data_path = candidates[0]
                name = entry
            else:
                data_path = path
                name = strip_typed_suffix(entry)

            data = read_array(data_path, default_dtype=idx_t).astype(idx_t)
            nodesets.append({"name": name, "nodes": np.unique(data)})

        return nodesets

    for sideset in sidesets:
        if not sideset["nodes"]:
            continue
        combined = np.concatenate(sideset["nodes"])
        nodesets.append({"name": sideset["name"], "nodes": np.unique(combined)})

    return nodesets


def load_time_whole(folder):
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


def write_field_names(dataset, var_name, dim_name, names):
    if not names:
        return None
    dataset.createDimension(dim_name, len(names))
    var = dataset.createVariable(var_name, "S1", (dim_name, "len_string"))
    var[:, :] = string_matrix(names)
    return var


def write_exodus(
    output_mesh,
    title,
    points,
    connectivity,
    element_type,
    blocks,
    sidesets,
    nodesets,
    point_fields,
    cell_fields,
    time_whole,
):
    import netCDF4

    n_dim, n_nodes = points.shape
    n_elem, nnodes_per_elem = connectivity.shape

    _, element_info = normalize_element_type(element_type)
    if nnodes_per_elem != element_info["nnodes"]:
        raise RuntimeError(
            f"connectivity width {nnodes_per_elem} does not match {element_type} ({element_info['nnodes']})"
        )

    time_whole = np.asarray(time_whole, dtype=np.float32)
    n_time_steps = len(time_whole)

    point_field_names = list(point_fields.keys())
    cell_field_names = list(cell_fields.keys())

    with netCDF4.Dataset(output_mesh, "w", format="NETCDF3_64BIT_OFFSET") as nc:
        nc.title = title
        nc.version = np.float32(5.1)
        nc.api_version = np.float32(5.1)
        nc.floating_point_word_size = np.int64(8)

        nc.createDimension("num_dim", n_dim)
        nc.createDimension("num_nodes", n_nodes)
        nc.createDimension("num_elem", n_elem)
        nc.createDimension("num_el_blk", len(blocks))
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

        coord = nc.createVariable("coord", points.dtype, ("num_dim", "num_nodes"))
        coord[:, :] = points

        eb_status = nc.createVariable("eb_status", np.int32, ("num_el_blk",))
        eb_status[:] = np.ones(len(blocks), dtype=np.int32)

        eb_prop1 = nc.createVariable("eb_prop1", np.int32, ("num_el_blk",))
        eb_prop1.setncattr("name", "ID")
        eb_prop1[:] = np.arange(1, len(blocks) + 1, dtype=np.int32)

        eb_names = nc.createVariable("eb_names", "S1", ("num_el_blk", "len_string"))
        eb_names[:, :] = string_matrix([block["name"] for block in blocks])

        for block_index, block in enumerate(blocks, start=1):
            block_size = block["end"] - block["begin"]
            nc.createDimension(f"num_el_in_blk{block_index}", block_size)
            nc.createDimension(f"num_nod_per_el{block_index}", nnodes_per_elem)
            connect = nc.createVariable(
                f"connect{block_index}",
                np.int32,
                (f"num_el_in_blk{block_index}", f"num_nod_per_el{block_index}"),
            )
            connect.setncattr("elem_type", element_info["exodus"])
            block_connect = connectivity[block["begin"] : block["end"], :] + 1
            connect[:, :] = block_connect.astype(np.int32)

        if len(blocks) > 1:
            for block_index, block in enumerate(blocks, start=1):
                prop = nc.createVariable(f"eb_prop{block_index + 1}", np.int32, ("num_el_blk",))
                prop.setncattr("name", block["name"])
                values = np.zeros(len(blocks), dtype=np.int32)
                values[block_index - 1] = block_index
                prop[:] = values

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
            elem_var_tab[:, :] = np.ones((len(blocks), len(cell_field_names)), dtype=np.int32)

            for var_index, name in enumerate(cell_field_names, start=1):
                first = cell_fields[name][0]
                for block_index, block in enumerate(blocks, start=1):
                    values = nc.createVariable(
                        f"vals_elem_var{var_index}eb{block_index}",
                        first.dtype,
                        ("time_step", f"num_el_in_blk{block_index}"),
                    )
                    for time_index in range(n_time_steps):
                        values[time_index, :] = cell_fields[name][time_index][
                            block["begin"] : block["end"]
                        ]

        if sidesets:
            ss_status = nc.createVariable("ss_status", np.int32, ("num_side_sets",))
            ss_status[:] = np.ones(len(sidesets), dtype=np.int32)

            ss_prop1 = nc.createVariable("ss_prop1", np.int32, ("num_side_sets",))
            ss_prop1.setncattr("name", "ID")
            ss_prop1[:] = np.arange(1, len(sidesets) + 1, dtype=np.int32)

            ss_names = nc.createVariable("ss_names", "S1", ("num_side_sets", "len_string"))
            ss_names[:, :] = string_matrix([sideset["name"] for sideset in sidesets])

            for ss_index, sideset in enumerate(sidesets, start=1):
                size = len(sideset["parent"])
                if len(sideset["lfi"]) != size:
                    raise RuntimeError(f"sideset '{sideset['name']}' has inconsistent lengths")

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
                elem_ss[:] = sideset["parent"].astype(np.int64) + 1
                side_ss[:] = sideset["lfi"].astype(np.int16) + 1

        if nodesets:
            ns_status = nc.createVariable("ns_status", np.int32, ("num_node_sets",))
            ns_status[:] = np.ones(len(nodesets), dtype=np.int32)

            ns_prop1 = nc.createVariable("ns_prop1", np.int32, ("num_node_sets",))
            ns_prop1.setncattr("name", "ID")
            ns_prop1[:] = np.arange(1, len(nodesets) + 1, dtype=np.int32)

            ns_names = nc.createVariable("ns_names", "S1", ("num_node_sets", "len_string"))
            ns_names[:, :] = string_matrix([nodeset["name"] for nodeset in nodesets])

            for ns_index, nodeset in enumerate(nodesets, start=1):
                nodes = np.unique(nodeset["nodes"].astype(np.int64))
                nc.createDimension(f"num_nod_ns{ns_index}", len(nodes))
                node_ns = nc.createVariable(
                    f"node_ns{ns_index}",
                    np.int32,
                    (f"num_nod_ns{ns_index}",),
                )
                node_ns[:] = nodes + 1


def raw_to_exodusII(input_folder, output_mesh, title=None):
    meta = read_simple_meta(os.path.join(input_folder, "meta.yaml"))
    element_type = meta.get("element_type")
    points = load_points(input_folder)
    connectivity = load_connectivity(input_folder)
    blocks = load_block_ranges(input_folder, connectivity.shape[0])
    sidesets = load_sidesets(input_folder)
    nodesets = load_nodesets(input_folder, sidesets)

    point_fields = load_field_series(os.path.join(input_folder, "point_data"), points.shape[1])
    cell_fields = load_field_series(os.path.join(input_folder, "cell_data"), connectivity.shape[0])
    time_whole = load_time_whole(input_folder)
    time_whole, _ = align_time_series(point_fields, cell_fields, time_whole)

    if title is None:
        title = f"Created by raw_to_exodusII.py from {os.path.basename(os.path.abspath(input_folder))}"

    write_exodus(
        output_mesh=output_mesh,
        title=title,
        points=points,
        connectivity=connectivity,
        element_type=element_type,
        blocks=blocks,
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
