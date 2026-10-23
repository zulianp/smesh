#!/usr/bin/env python3
"""Map nodal fields computed on part of a mesh back onto the whole mesh.

usage: raw_point_data_to_parent.py <fields> <part_mesh> <parent_mesh> <output> [--step=N | --all-steps]
                                   [--fields=a,b]

  <fields>       folder of field files: `<name>.<component>.<step>.<dtype>` (vectors, one file per
                 component) or `<name>.<step>.<dtype>` (scalars), as SFEM's Output writes them
  <part_mesh>    the mesh the fields live on; its `node_mapping.<dtype>` gives the parent node of
                 each of its nodes (a block-selective read writes it). Without one the part is the
                 parent itself.
  <parent_mesh>  the whole mesh, e.g. the folder exodusII_to_raw.py wrote
  <output>       a raw mesh folder: the parent's files (linked) and point_data/ with every field
                 on all parent nodes, zero on the nodes the part does not have. Vector components
                 are named <name>_x, _y, _z, which ParaView reads as one vector.

  --step=N       the step to map; the last one by default
  --all-steps    every step, as a time series: <name>_x.<k> for step k, and time_whole from the
                 fields folder's time.txt, which raw_to_exodusII.py writes as the file's times
  --fields=a,b   only these fields

Then `raw_to_exodusII.py <output> out.e` writes the Exodus file.
"""

import getopt
import os
import sys

import numpy as np

from common.raw_io import (
    find_stem_file,
    load_points,
    mkdir,
    read_array,
    select_steps,
    step_field_components,
    step_field_files,
    write_array,
)


def field_files(folder, wanted):
    """{name: {step: {component: path}}}, component None for scalars."""
    found, _ = step_field_files(os.path.join(folder, entry) for entry in os.listdir(folder))
    return {name: steps for name, steps in found.items() if not wanted or name in wanted}


def link_parent(parent, output):
    mkdir(output)
    for entry in os.listdir(parent):
        if entry in ("point_data", "cell_data"):
            continue
        target = os.path.join(output, entry)
        if os.path.lexists(target):
            os.remove(target)
        os.symlink(os.path.abspath(os.path.join(parent, entry)), target)


def main(argv):
    if len(argv) < 5:
        print(__doc__)
        sys.exit(1)
    fields_folder, part, parent, output = argv[1:5]
    step = None
    all_steps = False
    wanted = None
    try:
        opts, _ = getopt.getopt(argv[5:], "", ["step=", "all-steps", "fields="])
    except getopt.GetoptError as err:
        print(err)
        sys.exit(1)
    for opt, arg in opts:
        if opt == "--step":
            step = int(arg)
        elif opt == "--all-steps":
            all_steps = True
        elif opt == "--fields":
            wanted = set(arg.split(","))

    n_parent = load_points(parent)[0].shape[1]
    mapping_path = find_stem_file(part, "node_mapping")
    if mapping_path is not None:
        mapping = read_array(mapping_path).astype(np.int64)
    else:
        mapping = np.arange(load_points(part)[0].shape[1], dtype=np.int64)
    n_part = len(mapping)
    if n_part > n_parent or (mapping_path is None and n_part != n_parent):
        sys.exit(f"raw_point_data_to_parent: {part} has {n_part} nodes, the parent {n_parent}")
    if n_part and (mapping.min() < 0 or mapping.max() >= n_parent or len(np.unique(mapping)) != n_part):
        sys.exit(f"raw_point_data_to_parent: {mapping_path} is not a map into the {n_parent} parent nodes")

    found = field_files(fields_folder, wanted)
    if not found:
        sys.exit(f"raw_point_data_to_parent: no field files in {fields_folder}")
    if wanted and set(found) != wanted:
        sys.exit(f"raw_point_data_to_parent: missing fields {sorted(wanted - set(found))}")

    link_parent(parent, output)
    point_data = os.path.join(output, "point_data")
    mkdir(point_data)
    series = None
    for name, steps in sorted(found.items()):
        try:
            selected = select_steps(name, steps, step, all_steps)
            if all_steps and selected != list(range(len(selected))):
                raise RuntimeError(f"{name} has steps {selected[0]} to {selected[-1]} with gaps")
            if series is not None and all_steps and selected != series:
                raise RuntimeError(f"{name} has other steps than the fields before it")
            series = selected
            for s in selected:
                components = step_field_components(name, s, steps[s])
                # A time series names each step; a single step keeps the plain name.
                tag = f".{s}" if all_steps else ""
                for suffix, path in components:
                    values = read_array(path)
                    if len(values) != n_part:
                        sys.exit(f"raw_point_data_to_parent: {path} has {len(values)} values for {n_part} nodes")
                    whole = np.zeros(n_parent, dtype=values.dtype)
                    whole[mapping] = values
                    write_array(os.path.join(point_data, name + suffix + tag), whole, values.dtype)
        except RuntimeError as err:
            sys.exit(f"raw_point_data_to_parent: {err}")
        what = f"steps 0 to {selected[-1]}" if all_steps else f"step {selected[0]}"
        print(f"{name}: {what}, {len(components)} component(s), {n_part} of {n_parent} nodes")
    if all_steps:
        times_path = os.path.join(fields_folder, "time.txt")
        if not os.path.isfile(times_path):
            sys.exit(f"raw_point_data_to_parent: --all-steps needs {times_path}")
        times = np.loadtxt(times_path, dtype=np.float64).reshape(-1)
        if len(times) < len(series):
            sys.exit(f"raw_point_data_to_parent: {times_path} has {len(times)} times for {len(series)} steps")
        write_array(os.path.join(output, "time_whole"), times[: len(series)], np.float64)


if __name__ == "__main__":
    main(sys.argv)
