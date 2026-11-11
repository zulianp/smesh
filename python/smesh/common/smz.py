"""Compressed smesh arrays (`<stem>.<dtype>.smz`, see docs/compression.md).

read_info(path)   header fields (pure Python, no dependency on the C++ tools)
read_array(path)  the values as a numpy array; decoding runs the `smz` tool built with smesh
                  (taken from $SMESH_SMZ, or found on PATH): `smz decompress` into a temporary
                  raw file, read with numpy.
load(path, dtype) np.fromfile for raw files, read_array for .smz files
"""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
import tempfile

import numpy as np

MAGIC = b"\x89SMZ\r\n\x1a\n"
_DTYPES = {1: "int8", 2: "uint8", 3: "int16", 4: "uint16", 5: "int32", 6: "uint32", 7: "int64", 8: "uint64",
           9: "float16", 10: "float32", 11: "float64"}
_MODES = {0: "lossless", 1: "abs", 2: "range", 3: "pwrel", 4: "h"}


def is_smz(path):
    return str(path).endswith(".smz")


def typed_extension(path):
    """'float32' for x.float32 and x.float32.smz."""
    name = os.path.basename(str(path))
    if name.endswith(".smz"):
        name = name[: -len(".smz")]
    return name.rsplit(".", 1)[-1] if "." in name else ""


def read_info(path):
    with open(path, "rb") as f:
        h = f.read(128)
    if len(h) < 128 or h[:8] != MAGIC:
        raise RuntimeError(f"{path}: not an .smz file")
    major, minor = struct.unpack_from("<HH", h, 8)
    if major != 1:
        raise RuntimeError(f"{path}: unsupported .smz version {major}.{minor}")
    dtype, mode = h[16], h[17]
    n_values, chunk_values, n_chunks, payload = struct.unpack_from("<QQQQ", h, 24)
    eb, eb_abs, origin = struct.unpack_from("<ddd", h, 56)
    vmin, vmax, max_err = struct.unpack_from("<ddd", h, 88)
    (max_ratio,) = struct.unpack_from("<f", h, 112)
    file_bytes = 128 + 32 * n_chunks + payload
    return {
        "version": f"{major}.{minor}",
        "dtype": np.dtype(_DTYPES[dtype]),
        "mode": _MODES.get(mode, "?"),
        "eb": eb,
        "eb_abs": eb_abs,
        "n_values": n_values,
        "chunk_values": chunk_values,
        "n_chunks": n_chunks,
        "file_bytes": file_bytes,
        "ratio": n_values * np.dtype(_DTYPES[dtype]).itemsize / file_bytes if file_bytes else 0.0,
        "min": vmin,
        "max": vmax,
        "max_abs_err": max_err,
        "max_err_over_bound": max_ratio,
    }


def _tool():
    tool = os.environ.get("SMESH_SMZ") or shutil.which("smz")
    if not tool:
        raise RuntimeError("reading .smz files needs the `smz` tool built with smesh: put it on PATH or set SMESH_SMZ")
    return tool


def read_array(path):
    info = read_info(path)
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "values." + info["dtype"].name)
        res = subprocess.run([_tool(), "decompress", str(path), out], capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"smz decompress {path} failed: {res.stderr.strip()}")
        data = np.fromfile(out, dtype=info["dtype"])
    if len(data) != info["n_values"]:
        raise RuntimeError(f"{path}: expected {info['n_values']} values, got {len(data)}")
    return data


def load(path, dtype=None):
    """np.fromfile for raw files; .smz files are decoded (and converted to dtype if given)."""
    if is_smz(path):
        data = read_array(path)
        return data if dtype is None else data.astype(dtype, copy=False)
    return np.fromfile(path, dtype=dtype)
