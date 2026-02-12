import glob

def extension_to_dtype(extension):
    import numpy as np
    if extension == "raw":
        return np.float32
    return np.dtype(extension)
    
def extension(path):
    return path.split('.')[-1]

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