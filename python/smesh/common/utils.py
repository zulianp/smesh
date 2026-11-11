import glob
import numpy as np

def dtype_to_extension(dtype):
    return f"{np.dtype(dtype).name}"

def extension_to_dtype(extension):
    if extension == "raw":
        return np.float32
    return np.dtype(extension)

def extension(path):
    """The type extension: 'float32' for x.float32 and for compressed x.float32.smz."""
    parts = path.split('.')
    if len(parts) >= 3 and parts[-1] == 'smz':
        return parts[-2]
    return parts[-1]

def detect_files(pattern, extensions):
    files = glob.glob(pattern)
    paths = []
    for file in files:
        ext = extension(file)
        if ext in extensions:
            paths.append(file)
    return paths

def index_from_filename(filename):
    return int(filename.split('.')[0].split('i')[-1])
