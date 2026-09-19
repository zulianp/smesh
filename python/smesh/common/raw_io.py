"""C++ Mesh::write/read layout for Python converters.

Single-block: root i*.<idx>, x/y/z.<geom>, meta.yaml
Multi-block:  blocks/<name>/i*.<idx>, shared coords, meta.yaml with blocks:
Sidesets:     sidesets/<name>/{parent,lfi,meta.yaml}
              or sidesets/<name>/<block_id>/ when a name spans blocks
Nodesets:     nodesets/<name>/{nodes,meta.yaml}

parent is 0-based local element index in block_id; lfi is 0-based (Exodus side - 1).
"""

from __future__ import annotations

import os
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from common.utils import detect_files, dtype_to_extension

EXODUS_OUTPUT_EXTENSIONS = (".e", ".exo", ".ex2")
INDEX_EXTENSIONS = ("raw", "int16", "int32", "int64")
GEOM_EXTENSIONS = ("raw", "float16", "float32", "float64")
DEFAULT_IDX_DTYPE = np.dtype(np.int32)
DEFAULT_GEOM_DTYPE = np.dtype(np.float64)
DEFAULT_PARENT_DTYPE = np.dtype(np.int32)
DEFAULT_LFI_DTYPE = np.dtype(np.int16)
DEFAULT_GEOM_MAP = "IsoParametric"
MAX_NODES_X_ELEMENT = 27

# C++ type_to_string names -> Exodus elem_type attribute and node count.
ELEMENT_INFO = {
    "NODE1": {"smesh": "NODE1", "exodus": "NODE", "nnodes": 1},
    "EDGE2": {"smesh": "EDGE2", "exodus": "BAR2", "nnodes": 2},
    "EDGE3": {"smesh": "EDGE3", "exodus": "BAR3", "nnodes": 3},
    "TRI3": {"smesh": "TRI3", "exodus": "TRI3", "nnodes": 3},
    "TRISHELL3": {"smesh": "TRISHELL3", "exodus": "TRI3", "nnodes": 3},
    "TRI6": {"smesh": "TRI6", "exodus": "TRI6", "nnodes": 6},
    "TRISHELL6": {"smesh": "TRISHELL6", "exodus": "TRI6", "nnodes": 6},
    "QUAD4": {"smesh": "QUAD4", "exodus": "QUAD4", "nnodes": 4},
    "QUADSHELL4": {"smesh": "QUADSHELL4", "exodus": "SHELL4", "nnodes": 4},
    "QUAD9": {"smesh": "QUAD9", "exodus": "QUAD9", "nnodes": 9},
    "QUADSHELL9": {"smesh": "QUADSHELL9", "exodus": "SHELL9", "nnodes": 9},
    "TET4": {"smesh": "TET4", "exodus": "TETRA", "nnodes": 4},
    "TET10": {"smesh": "TET10", "exodus": "TETRA10", "nnodes": 10},
    "HEX8": {"smesh": "HEX8", "exodus": "HEX8", "nnodes": 8},
    "HEX27": {"smesh": "HEX27", "exodus": "HEX27", "nnodes": 27},
    "PROTEUS_HEX27": {"smesh": "PROTEUS_HEX27", "exodus": "HEX27", "nnodes": 27},
    "WEDGE6": {"smesh": "WEDGE6", "exodus": "WEDGE", "nnodes": 6},
    "PYRAMID5": {"smesh": "PYRAMID5", "exodus": "PYRAMID", "nnodes": 5},
}

ELEMENT_ALIASES = {
    "HEX": "HEX8",
    "HEX8": "HEX8",
    "hex": "HEX8",
    "hex8": "HEX8",
    "hexahedron": "HEX8",
    "hexahedron8": "HEX8",
    "HEX27": "HEX27",
    "hex27": "HEX27",
    "hexahedron27": "HEX27",
    "PROTEUS_HEX27": "PROTEUS_HEX27",
    "proteus_hex27": "PROTEUS_HEX27",
    "TET4": "TET4",
    "TETRA": "TET4",
    "TETRA4": "TET4",
    "tetra": "TET4",
    "tetra4": "TET4",
    "TET10": "TET10",
    "TETRA10": "TET10",
    "tetra10": "TET10",
    "QUAD4": "QUAD4",
    "QUAD": "QUAD4",
    "quad": "QUAD4",
    "quad4": "QUAD4",
    "SHELL": "QUADSHELL4",
    "SHELL4": "QUADSHELL4",
    "QUADSHELL4": "QUADSHELL4",
    "quadshell4": "QUADSHELL4",
    "QUAD9": "QUAD9",
    "quad9": "QUAD9",
    "SHELL9": "QUADSHELL9",
    "QUADSHELL9": "QUADSHELL9",
    "TRI3": "TRI3",
    "TRI": "TRI3",
    "TRIANGLE": "TRI3",
    "tri": "TRI3",
    "tri3": "TRI3",
    "triangle": "TRI3",
    "TRISHELL3": "TRISHELL3",
    "TRI6": "TRI6",
    "tri6": "TRI6",
    "triangle6": "TRI6",
    "TRISHELL6": "TRISHELL6",
    "WEDGE6": "WEDGE6",
    "WEDGE": "WEDGE6",
    "wedge": "WEDGE6",
    "wedge6": "WEDGE6",
    "prism": "WEDGE6",
    "prism6": "WEDGE6",
    "PENTA": "WEDGE6",
    "PYRAMID5": "PYRAMID5",
    "PYRAMID": "PYRAMID5",
    "pyramid": "PYRAMID5",
    "pyramid5": "PYRAMID5",
}

MESHIO_TO_SMESH = {
    "vertex": "NODE1",
    "line": "EDGE2",
    "line3": "EDGE3",
    "triangle": "TRI3",
    "triangle6": "TRI6",
    "quad": "QUAD4",
    "quad9": "QUAD9",
    "tetra": "TET4",
    "tetra10": "TET10",
    "hexahedron": "HEX8",
    "hexahedron27": "HEX27",
    "wedge": "WEDGE6",
    "pyramid": "PYRAMID5",
}

SMESH_TO_MESHIO = {
    "NODE1": "vertex",
    "EDGE2": "line",
    "EDGE3": "line3",
    "TRI3": "triangle",
    "TRISHELL3": "triangle",
    "TRI6": "triangle6",
    "TRISHELL6": "triangle6",
    "QUAD4": "quad",
    "QUADSHELL4": "quad",
    "QUAD9": "quad9",
    "QUADSHELL9": "quad9",
    "TET4": "tetra",
    "TET10": "tetra10",
    "HEX8": "hexahedron",
    "HEX27": "hexahedron27",
    "PROTEUS_HEX27": "hexahedron27",
    "WEDGE6": "wedge",
    "PYRAMID5": "pyramid",
}

ELEMENT_TYPE_BY_NUM_NODES = {
    1: "NODE1",
    2: "EDGE2",
    3: "TRI3",
    4: "TET4",
    5: "PYRAMID5",
    6: "WEDGE6",
    8: "HEX8",
    9: "QUAD9",
    10: "TET10",
    27: "HEX27",
}


def is_exodus_path(path):
    return os.path.splitext(path)[1].lower() in EXODUS_OUTPUT_EXTENSIONS


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
    if default_dtype is not None:
        return np.dtype(default_dtype)
    return None


def read_array(path, default_dtype=None, count=None):
    dtype = dtype_from_path(path, default_dtype)
    if dtype is None:
        if count is not None:
            nbytes = os.path.getsize(path)
            if count == 0:
                dtype = np.dtype(default_dtype if default_dtype is not None else np.int32)
            else:
                itemsize = nbytes // count if count else 0
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


def write_array(path_stem, data, dtype):
    dtype = np.dtype(dtype)
    path = f"{path_stem}.{dtype_to_extension(dtype)}"
    np.ascontiguousarray(np.asarray(data, dtype=dtype)).tofile(path)
    return path


def mkdir(path):
    if path and not os.path.exists(path):
        os.makedirs(path)


def find_stem_file(folder, stem):
    if not os.path.isdir(folder):
        return None
    matches = []
    for entry in os.listdir(folder):
        path = os.path.join(folder, entry)
        if not os.path.isfile(path):
            continue
        if strip_typed_suffix(entry) == stem:
            matches.append(path)
    matches.sort()
    return matches[0] if matches else None


def smesh_element_type(raw_type, nnodes=None, spatial_dim=3):
    if raw_type is not None:
        key = str(raw_type).strip()
        if key in ELEMENT_ALIASES:
            return ELEMENT_ALIASES[key]
        upper = key.upper()
        if upper in ELEMENT_ALIASES:
            return ELEMENT_ALIASES[upper]
        if key in ELEMENT_INFO:
            return key
        raise RuntimeError(f"unsupported element_type '{raw_type}'")
    if nnodes is None:
        raise RuntimeError("missing element_type and connectivity width")
    nnodes = int(nnodes)
    if nnodes == 4:
        return "QUAD4" if spatial_dim == 2 else "TET4"
    if nnodes == 6:
        return "TRI6" if spatial_dim == 2 else "WEDGE6"
    inferred = ELEMENT_TYPE_BY_NUM_NODES.get(nnodes)
    if inferred is None:
        raise RuntimeError(f"unable to infer element_type from connectivity width {nnodes}")
    return inferred


def exodus_element_type(smesh_type):
    info = ELEMENT_INFO.get(smesh_element_type(smesh_type))
    if info is None:
        raise RuntimeError(f"unsupported element_type '{smesh_type}'")
    return info["exodus"]


def element_nnodes(smesh_type):
    info = ELEMENT_INFO.get(smesh_element_type(smesh_type))
    if info is None:
        raise RuntimeError(f"unsupported element_type '{smesh_type}'")
    return info["nnodes"]


def meshio_cell_type(smesh_type):
    mapped = smesh_element_type(smesh_type)
    if mapped not in SMESH_TO_MESHIO:
        raise RuntimeError(f"no meshio cell type for '{smesh_type}'")
    return SMESH_TO_MESHIO[mapped]


def smesh_type_from_meshio(cell_type, nnodes=None, spatial_dim=3):
    if cell_type in MESHIO_TO_SMESH:
        name = MESHIO_TO_SMESH[cell_type]
        if cell_type == "quad" and nnodes == 9:
            return "QUAD9"
        if cell_type == "triangle" and nnodes == 6:
            return "TRI6"
        if cell_type == "tetra" and nnodes == 10:
            return "TET10"
        if cell_type == "hexahedron" and nnodes == 27:
            return "HEX27"
        if cell_type == "quad" and spatial_dim == 3 and nnodes == 4:
            return "QUAD4"
        return name
    return smesh_element_type(cell_type, nnodes=nnodes)


def unique_preserve_order(nodes):
    nodes = np.asarray(nodes)
    if nodes.size == 0:
        return nodes
    _, first = np.unique(nodes, return_index=True)
    if first.size == nodes.size:
        return nodes
    return nodes[np.sort(first)]


@dataclass
class Block:
    name: str
    element_type: str
    connectivity: np.ndarray
    geom_map: str = DEFAULT_GEOM_MAP
    exodus_id: Optional[int] = None

    @property
    def n_elements(self):
        if self.connectivity.ndim != 2:
            return 0
        return int(self.connectivity.shape[0])

    @property
    def nnodes(self):
        if self.connectivity.ndim != 2:
            return 0
        return int(self.connectivity.shape[1])


@dataclass
class Sideset:
    name: str
    block_id: int
    parent: np.ndarray
    lfi: np.ndarray
    exodus_id: Optional[int] = None

    @property
    def size(self):
        return int(len(self.parent))


@dataclass
class Nodeset:
    name: str
    nodes: np.ndarray
    exodus_id: Optional[int] = None

    @property
    def size(self):
        return int(len(self.nodes))


@dataclass
class RawMesh:
    points: np.ndarray
    blocks: List[Block] = field(default_factory=list)
    sidesets: List[Sideset] = field(default_factory=list)
    nodesets: List[Nodeset] = field(default_factory=list)
    idx_dtype: np.dtype = DEFAULT_IDX_DTYPE
    geom_dtype: np.dtype = DEFAULT_GEOM_DTYPE
    parent_dtype: np.dtype = DEFAULT_PARENT_DTYPE

    @property
    def n_nodes(self):
        if self.points.ndim != 2:
            return 0
        return int(self.points.shape[1])

    @property
    def spatial_dimension(self):
        if self.points.ndim != 2:
            return 0
        return int(self.points.shape[0])

    @property
    def n_elements(self):
        return int(sum(block.n_elements for block in self.blocks))

    def element_offsets(self):
        nels = [block.n_elements for block in self.blocks]
        return np.cumsum(np.asarray([0] + nels, dtype=np.int64))


def element_offsets_from_counts(counts):
    return np.cumsum(np.asarray([0] + list(counts), dtype=np.int64))


def global_to_block_local(global_0based, offsets):
    """Map 0-based global element ids to (block_id, local) via prefix sums."""
    elem = np.asarray(global_0based, dtype=np.int64)
    if elem.size == 0:
        return (
            np.zeros((0,), dtype=np.int32),
            np.zeros((0,), dtype=np.int64),
        )
    block_id = np.searchsorted(offsets[1:], elem, side="right").astype(np.int32)
    local = elem - offsets[block_id]
    return block_id, local


def block_local_to_global(block_id, local, offsets):
    block_id = np.asarray(block_id, dtype=np.int64)
    local = np.asarray(local, dtype=np.int64)
    return offsets[block_id] + local


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
            if key in ("name", "element_type", "cell_type", "elem_type", "geom_map"):
                current[key] = value
            elif key in ("elem_num_nodes", "n_elements", "exodus_id", "block_id"):
                current[key] = int(value)
    if current:
        blocks.append(current)
    return blocks


def parse_points_meta(meta_path):
    """Return {axis: filename} from meta.yaml `points:` if present."""
    names = {}
    if not os.path.exists(meta_path):
        return names
    in_points = False
    with open(meta_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if stripped == "points:":
                in_points = True
                continue
            if not in_points:
                continue
            if (
                not line.startswith(" ")
                and not line.startswith("\t")
                and not stripped.startswith("-")
            ):
                break
            if stripped.startswith("- ") and ":" in stripped:
                item = stripped[2:]
                key, value = item.split(":", 1)
                key = key.strip()
                value = value.strip()
                if key in ("x", "y", "z") and value:
                    names[key] = value
    return names


_GEOM_EXT_RANK = {"float64": 0, "float32": 1, "float16": 2, "raw": 3}


def _pick_axis_file(folder, axis, listed=None):
    if listed:
        path = listed if os.path.isabs(listed) else os.path.join(folder, listed)
        if os.path.isfile(path):
            return path
    matches = []
    if not os.path.isdir(folder):
        return None
    for entry in os.listdir(folder):
        path = os.path.join(folder, entry)
        if os.path.isfile(path) and strip_typed_suffix(entry) == axis:
            matches.append(path)
    matches.sort(
        key=lambda p: (
            _GEOM_EXT_RANK.get(os.path.basename(p).rsplit(".", 1)[-1], 9),
            p,
        )
    )
    return matches[0] if matches else None


def _remove_stale_stem(folder, stem, keep_path):
    if not os.path.isdir(folder):
        return
    keep = os.path.abspath(keep_path) if keep_path else None
    for entry in os.listdir(folder):
        path = os.path.join(folder, entry)
        if not os.path.isfile(path):
            continue
        if strip_typed_suffix(entry) != stem:
            continue
        if keep is not None and os.path.abspath(path) == keep:
            continue
        os.remove(path)


def load_points(folder):
    listed = parse_points_meta(os.path.join(folder, "meta.yaml"))
    points = []
    geom_dtype = None
    for axis in ("x", "y", "z"):
        path = _pick_axis_file(folder, axis, listed.get(axis))
        if path is None:
            path_list = detect_files(f"{folder}/{axis}.*", list(GEOM_EXTENSIONS))
            if path_list:
                path_list.sort(
                    key=lambda p: (
                        _GEOM_EXT_RANK.get(os.path.basename(p).rsplit(".", 1)[-1], 9),
                        p,
                    )
                )
                path = path_list[0]
        if path is None:
            break
        dtype = dtype_from_path(path, DEFAULT_GEOM_DTYPE)
        axis_data = np.fromfile(path, dtype=dtype)
        points.append(axis_data)
        geom_dtype = dtype
    if not points:
        raise RuntimeError(f"no coordinate files found in {folder}")
    n_nodes = len(points[0])
    for axis in points:
        if len(axis) != n_nodes:
            raise RuntimeError("coordinate arrays have inconsistent lengths")
    return np.vstack(points), np.dtype(geom_dtype)


def write_points(folder, points, geom_dtype):
    mkdir(folder)
    points = np.asarray(points)
    geom_dtype = np.dtype(geom_dtype)
    axes = ("x", "y", "z")
    for d in range(points.shape[0]):
        keep = write_array(os.path.join(folder, axes[d]), points[d, :], geom_dtype)
        _remove_stale_stem(folder, axes[d], keep)


def load_connectivity(folder):
    index_pattern = re.compile(r"^i(\d+)\.")
    by_index = {}
    if not os.path.isdir(folder):
        raise RuntimeError(f"no connectivity folder {folder}")
    for entry in os.listdir(folder):
        match = index_pattern.match(entry)
        if match:
            by_index.setdefault(int(match.group(1)), []).append(os.path.join(folder, entry))
    if not by_index:
        raise RuntimeError(f"no connectivity files found in {folder}")
    preferred_ext = dtype_to_extension(DEFAULT_IDX_DTYPE)
    arrays = []
    idx_dtype = DEFAULT_IDX_DTYPE
    for ii in range(max(by_index) + 1):
        paths = by_index.get(ii)
        if not paths:
            raise RuntimeError(f"missing connectivity file i{ii} in {folder}")
        path = None
        for candidate in paths:
            if candidate.endswith("." + preferred_ext):
                path = candidate
                break
        if path is None:
            path = sorted(paths)[0]
        dtype = dtype_from_path(path, DEFAULT_IDX_DTYPE)
        arrays.append(np.fromfile(path, dtype=dtype))
        idx_dtype = dtype
    n_elements = len(arrays[0])
    for array in arrays:
        if len(array) != n_elements:
            raise RuntimeError("connectivity arrays have inconsistent lengths")
    conn = np.column_stack(arrays) if arrays else np.zeros((0, 0), dtype=idx_dtype)
    return conn, np.dtype(idx_dtype)


def write_connectivity(folder, connectivity, idx_dtype):
    mkdir(folder)
    connectivity = np.asarray(connectivity)
    idx_dtype = np.dtype(idx_dtype)
    if connectivity.ndim != 2:
        raise RuntimeError("connectivity must be (n_elements, nnodes)")
    nxe = connectivity.shape[1]
    for d in range(nxe):
        keep = write_array(os.path.join(folder, f"i{d}"), connectivity[:, d], idx_dtype)
        _remove_stale_stem(folder, f"i{d}", keep)
    index_pattern = re.compile(r"^i(\d+)\.")
    for entry in os.listdir(folder):
        match = index_pattern.match(entry)
        if match and int(match.group(1)) >= nxe:
            os.remove(os.path.join(folder, entry))


def _write_single_block_meta(folder, mesh, idx_ext, geom_ext):
    block = mesh.blocks[0]
    nxe = block.nnodes
    path = os.path.join(folder, "meta.yaml")
    lines = [
        "# smesh mesh meta file",
        f"spatial_dimension: {mesh.spatial_dimension}",
        f"elem_num_nodes: {nxe}",
        f"element_type: {block.element_type}",
        f"geom_map: {block.geom_map or DEFAULT_GEOM_MAP}",
        f"n_elements: {block.n_elements}",
        f"n_nodes: {mesh.n_nodes}",
    ]
    if block.exodus_id is not None:
        lines.append(f"exodus_id: {int(block.exodus_id)}")
    lines.append("elements:")
    for d in range(nxe):
        lines.append(f"- i{d}: i{d}.{idx_ext}")
    lines.append("points:")
    axes = ("x", "y", "z")
    for d in range(mesh.spatial_dimension):
        lines.append(f"- {axes[d]}: {axes[d]}.{geom_ext}")
    lines.append("rpath: true")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _write_multi_block_meta(folder, mesh, idx_ext, geom_ext):
    path = os.path.join(folder, "meta.yaml")
    lines = [
        "# SMESH mesh meta file (generated by smesh_write.cpp)",
        f"spatial_dimension: {mesh.spatial_dimension}",
        f"n_blocks: {len(mesh.blocks)}",
        "blocks:",
    ]
    for block in mesh.blocks:
        nxe = block.nnodes
        lines.append(f"- name: {block.name}")
        lines.append(f"  element_type: {block.element_type}")
        lines.append(f"  geom_map: {block.geom_map or DEFAULT_GEOM_MAP}")
        if block.exodus_id is not None:
            lines.append(f"  exodus_id: {int(block.exodus_id)}")
        lines.append(f"  elem_num_nodes: {nxe}")
        lines.append(f"  n_elements: {block.n_elements}")
        lines.append("  elements:")
        for d in range(nxe):
            lines.append(f"  - i{d}: blocks/{block.name}/i{d}.{idx_ext}")
    lines.append(f"n_nodes: {mesh.n_nodes}")
    lines.append("points:")
    axes = ("x", "y", "z")
    for d in range(mesh.spatial_dimension):
        lines.append(f"- {axes[d]}: {axes[d]}.{geom_ext}")
    lines.append("rpath: true")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def write_sideset_meta(folder, size, block_id, parent_name, lfi_name, exodus_id=None):
    mkdir(folder)
    path = os.path.join(folder, "meta.yaml")
    lines = [
        "# Automatically generated by smesh_Sideset.cpp",
        f"size: {int(size)}",
        f"block_id: {int(block_id)}",
        f"parent: {parent_name}",
        f"lfi: {lfi_name}",
        "rpath: true",
    ]
    if exodus_id is not None:
        lines.append(f"exodus_id: {int(exodus_id)}")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def write_nodeset_meta(folder, size, nodes_name, exodus_id=None):
    mkdir(folder)
    path = os.path.join(folder, "meta.yaml")
    lines = [
        "# Automatically generated by smesh_nodeset.cpp",
        f"size: {int(size)}",
        f"nodes: {nodes_name}",
        "rpath: true",
    ]
    if exodus_id is not None:
        lines.append(f"exodus_id: {int(exodus_id)}")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _load_one_sideset(folder, name, default_block_id=0):
    parent_path = find_stem_file(folder, "parent")
    lfi_path = find_stem_file(folder, "lfi")
    if parent_path is None or lfi_path is None:
        return None
    lfi = read_array(lfi_path, default_dtype=DEFAULT_LFI_DTYPE).astype(np.int16, copy=False)
    parent = read_array(parent_path, default_dtype=DEFAULT_PARENT_DTYPE, count=len(lfi))
    meta = read_simple_meta(os.path.join(folder, "meta.yaml"))
    block_id = default_block_id
    if "block_id" in meta and meta["block_id"] != "":
        block_id = int(meta["block_id"])
    exodus_id = int(meta["exodus_id"]) if meta.get("exodus_id") else None
    return Sideset(
        name=name,
        block_id=block_id,
        parent=np.asarray(parent, dtype=parent.dtype),
        lfi=lfi,
        exodus_id=exodus_id,
    )


def load_sidesets(folder):
    sidesets_dir = os.path.join(folder, "sidesets")
    if not os.path.isdir(sidesets_dir):
        return []

    sidesets = []
    for name in sorted(os.listdir(sidesets_dir)):
        ss_dir = os.path.join(sidesets_dir, name)
        if not os.path.isdir(ss_dir):
            continue
        direct = _load_one_sideset(ss_dir, name, default_block_id=0)
        if direct is not None:
            sidesets.append(direct)
            continue
        for child in sorted(os.listdir(ss_dir)):
            child_dir = os.path.join(ss_dir, child)
            if not os.path.isdir(child_dir):
                continue
            default_block_id = 0
            if child.isdigit() or (child.startswith("-") and child[1:].isdigit()):
                default_block_id = int(child)
            member = _load_one_sideset(child_dir, name, default_block_id=default_block_id)
            if member is not None:
                sidesets.append(member)
    return sidesets


def load_nodesets(folder):
    nodesets_dir = os.path.join(folder, "nodesets")
    if not os.path.isdir(nodesets_dir):
        return []

    nodesets = []
    for entry in sorted(os.listdir(nodesets_dir)):
        path = os.path.join(nodesets_dir, entry)
        if os.path.isdir(path):
            data_path = find_stem_file(path, "nodes")
            if data_path is None:
                data_path = find_stem_file(path, "nodeset")
            if data_path is None:
                data_path = find_stem_file(path, entry)
            if data_path is None:
                files = [
                    os.path.join(path, child)
                    for child in sorted(os.listdir(path))
                    if os.path.isfile(os.path.join(path, child))
                    and not child.endswith(".yaml")
                ]
                data_path = files[0] if files else None
            if data_path is None:
                continue
            name = entry
            meta = read_simple_meta(os.path.join(path, "meta.yaml"))
        elif os.path.isfile(path):
            data_path = path
            name = strip_typed_suffix(entry)
            meta = {}
        else:
            continue
        data = read_array(data_path, default_dtype=DEFAULT_IDX_DTYPE)
        exodus_id = int(meta["exodus_id"]) if meta.get("exodus_id") else None
        nodesets.append(Nodeset(name=name, nodes=np.asarray(data), exodus_id=exodus_id))
    return nodesets


def _discover_block_dirs(mesh_folder):
    blocks_root = os.path.join(mesh_folder, "blocks")
    if not os.path.isdir(blocks_root):
        return []
    names = []
    for name in sorted(os.listdir(blocks_root)):
        folder = os.path.join(blocks_root, name)
        if not os.path.isdir(folder):
            continue
        if detect_files(f"{folder}/i0.*", list(INDEX_EXTENSIONS)):
            names.append(name)
    return names


def _load_legacy_block_ranges(blocks_dir, n_elements):
    if not os.path.isdir(blocks_dir):
        return []
    blocks = []
    for entry in sorted(os.listdir(blocks_dir)):
        path = os.path.join(blocks_dir, entry)
        if not os.path.isfile(path):
            continue
        values = read_array(path, default_dtype=np.int64)
        if len(values) != 2:
            continue
        begin = int(values[0])
        end = int(values[1])
        if begin < 0 or end < begin or end > n_elements:
            raise RuntimeError(f"invalid block bounds [{begin}, {end}) in {path}")
        blocks.append((strip_typed_suffix(entry), begin, end))
    return blocks


def load_raw_mesh(folder):
    meta_path = os.path.join(folder, "meta.yaml")
    meta = read_simple_meta(meta_path)
    points, geom_dtype = load_points(folder)
    block_specs = parse_smesh_block_list(meta_path)
    dir_names = _discover_block_dirs(folder)

    blocks = []
    idx_dtype = DEFAULT_IDX_DTYPE
    parent_dtype = DEFAULT_PARENT_DTYPE

    if block_specs or dir_names:
        names = [spec["name"] for spec in block_specs] if block_specs else dir_names
        spec_by_name = {spec["name"]: spec for spec in block_specs}
        for name in names:
            block_folder = os.path.join(folder, "blocks", name)
            conn, idx_dtype = load_connectivity(block_folder)
            spec = spec_by_name.get(name, {})
            raw_type = spec.get("element_type") or spec.get("cell_type") or spec.get("elem_type")
            element_type = smesh_element_type(
                raw_type,
                nnodes=conn.shape[1] if conn.size else None,
                spatial_dim=points.shape[0],
            )
            geom_map = spec.get("geom_map") or DEFAULT_GEOM_MAP
            exodus_id = spec.get("exodus_id")
            blocks.append(
                Block(
                    name=name,
                    element_type=element_type,
                    connectivity=conn,
                    geom_map=geom_map,
                    exodus_id=int(exodus_id) if exodus_id is not None else None,
                )
            )
    else:
        root_has_conn = bool(detect_files(f"{folder}/i0.*", list(INDEX_EXTENSIONS)))
        if not root_has_conn:
            raise RuntimeError(f"no connectivity in {folder}")
        conn, idx_dtype = load_connectivity(folder)
        legacy = _load_legacy_block_ranges(os.path.join(folder, "blocks"), conn.shape[0])
        if legacy:
            for name, begin, end in legacy:
                piece = conn[begin:end, :]
                element_type = smesh_element_type(
                    None, nnodes=piece.shape[1], spatial_dim=points.shape[0]
                )
                blocks.append(
                    Block(
                        name=name,
                        element_type=element_type,
                        connectivity=piece,
                        geom_map=DEFAULT_GEOM_MAP,
                    )
                )
        else:
            raw_type = meta.get("element_type") or meta.get("cell_type") or meta.get("elem_type")
            element_type = smesh_element_type(
                raw_type,
                nnodes=conn.shape[1] if conn.size else None,
                spatial_dim=points.shape[0],
            )
            geom_map = meta.get("geom_map") or DEFAULT_GEOM_MAP
            exodus_id = int(meta["exodus_id"]) if meta.get("exodus_id") else None
            blocks.append(
                Block(
                    name="default",
                    element_type=element_type,
                    connectivity=conn,
                    geom_map=geom_map,
                    exodus_id=exodus_id,
                )
            )

    sidesets = load_sidesets(folder)
    nodesets = load_nodesets(folder)
    if sidesets:
        parent_dtype = np.dtype(sidesets[0].parent.dtype)

    return RawMesh(
        points=points,
        blocks=blocks,
        sidesets=sidesets,
        nodesets=nodesets,
        idx_dtype=idx_dtype,
        geom_dtype=geom_dtype,
        parent_dtype=parent_dtype,
    )


def write_sidesets(folder, sidesets, parent_dtype=DEFAULT_PARENT_DTYPE, lfi_dtype=DEFAULT_LFI_DTYPE):
    if not sidesets:
        return
    parent_dtype = np.dtype(parent_dtype)
    lfi_dtype = np.dtype(lfi_dtype)
    groups = OrderedDict()
    for ss in sidesets:
        groups.setdefault(ss.name, []).append(ss)

    root = os.path.join(folder, "sidesets")
    mkdir(root)
    for name, members in groups.items():
        name_dir = os.path.join(root, name)
        mkdir(name_dir)
        nested = len(members) > 1
        for ss in members:
            dest = os.path.join(name_dir, str(int(ss.block_id))) if nested else name_dir
            mkdir(dest)
            parent_path = write_array(
                os.path.join(dest, "parent"), ss.parent, parent_dtype
            )
            lfi_path = write_array(os.path.join(dest, "lfi"), ss.lfi, lfi_dtype)
            write_sideset_meta(
                dest,
                size=len(ss.parent),
                block_id=ss.block_id,
                parent_name=os.path.basename(parent_path),
                lfi_name=os.path.basename(lfi_path),
                exodus_id=ss.exodus_id,
            )


def write_nodesets(folder, nodesets, idx_dtype=DEFAULT_IDX_DTYPE):
    if not nodesets:
        return
    idx_dtype = np.dtype(idx_dtype)
    root = os.path.join(folder, "nodesets")
    mkdir(root)
    for ns in nodesets:
        dest = os.path.join(root, ns.name)
        mkdir(dest)
        nodes_path = write_array(os.path.join(dest, "nodes"), ns.nodes, idx_dtype)
        write_nodeset_meta(
            dest,
            size=len(ns.nodes),
            nodes_name=os.path.basename(nodes_path),
            exodus_id=ns.exodus_id,
        )


def write_raw_mesh(folder, mesh):
    mkdir(folder)
    idx_dtype = np.dtype(mesh.idx_dtype)
    geom_dtype = np.dtype(mesh.geom_dtype)
    parent_dtype = np.dtype(mesh.parent_dtype)
    idx_ext = dtype_to_extension(idx_dtype)
    geom_ext = dtype_to_extension(geom_dtype)

    write_points(folder, mesh.points, geom_dtype)

    named_single = (
        len(mesh.blocks) == 1
        and mesh.blocks[0].name not in ("", "default")
    )
    if len(mesh.blocks) == 1 and not named_single:
        write_connectivity(folder, mesh.blocks[0].connectivity, idx_dtype)
        _write_single_block_meta(folder, mesh, idx_ext, geom_ext)
    else:
        for block in mesh.blocks:
            write_connectivity(
                os.path.join(folder, "blocks", block.name),
                block.connectivity,
                idx_dtype,
            )
        _write_multi_block_meta(folder, mesh, idx_ext, geom_ext)
        index_pattern = re.compile(r"^i\d+\.")
        for entry in os.listdir(folder):
            if index_pattern.match(entry):
                os.remove(os.path.join(folder, entry))

    write_sidesets(folder, mesh.sidesets, parent_dtype=parent_dtype)
    write_nodesets(folder, mesh.nodesets, idx_dtype=idx_dtype)


def group_sidesets_by_name(sidesets):
    groups = OrderedDict()
    for ss in sidesets:
        groups.setdefault(ss.name, []).append(ss)
    return groups


def sidesets_to_global(sidesets, offsets):
    """Merge per-block sidesets of the same name into Exodus global (elem, side)."""
    out = []
    for name, members in group_sidesets_by_name(sidesets).items():
        members = sorted(members, key=lambda item: item.block_id)
        if not members:
            continue
        parents = []
        lfis = []
        for ss in members:
            local = np.asarray(ss.parent, dtype=np.int64)
            parents.append(offsets[ss.block_id] + local)
            lfis.append(np.asarray(ss.lfi, dtype=np.int16))
        parent = np.concatenate(parents) if parents else np.zeros((0,), dtype=np.int64)
        lfi = np.concatenate(lfis) if lfis else np.zeros((0,), dtype=np.int16)
        out.append(
            {
                "name": name,
                "parent": parent,
                "lfi": lfi,
                "exodus_id": members[0].exodus_id,
            }
        )
    return out


def split_global_sideset(name, elem_0based, side_0based, offsets, exodus_id=None):
    """Split an Exodus sideset into per-block Sideset objects."""
    elem = np.asarray(elem_0based, dtype=np.int64)
    side = np.asarray(side_0based, dtype=np.int16)
    if elem.size == 0:
        return [
            Sideset(
                name=name,
                block_id=0,
                parent=np.zeros((0,), dtype=DEFAULT_PARENT_DTYPE),
                lfi=np.zeros((0,), dtype=np.int16),
                exodus_id=exodus_id,
            )
        ]
    block_id, local = global_to_block_local(elem, offsets)
    members = []
    for bid in np.unique(block_id):
        mask = block_id == bid
        members.append(
            Sideset(
                name=name,
                block_id=int(bid),
                parent=local[mask].astype(DEFAULT_PARENT_DTYPE, copy=False),
                lfi=side[mask],
                exodus_id=exodus_id,
            )
        )
    return members
