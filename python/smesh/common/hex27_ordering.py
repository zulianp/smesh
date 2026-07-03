"""HEX27 connectivity permutations shared by raw_to_db and raw_to_exodusII."""

PROTEUS_HEX27_NAMES = frozenset({"PROTEUS_HEX27", "proteus_hex27"})

# PROTEUS_HEX27: Cartesian iz*9+iy*3+ix layout.
proteus_hex27_to_exodus_hex27 = (
    0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23,
    25, 21, 9, 11, 17, 15, 10, 14, 16, 12, 4, 22, 13,
)

# SFEM/Exodus and VTK agree on corners (0-19) and body node (26); face centers 20-25 differ.
# Exodus/PATRAN: 20=-y, 21=+x, 22=+y, 23=-x, 24=-z, 25=+z  -> .e / IOSS / FEM
# VTK tri-quadratic hex: 20=-x, 21=+x, 22=-y, 23=+y, 24=-z, 25=+z -> .vtu only
exodus_hex27_to_vtk_hex27 = (
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    23, 21, 20, 22, 24, 25, 26,
)


def prepare_exodus_hex27_connectivity(connectivity, element_type):
    """PROTEUS_HEX27 -> HEX27; otherwise return connectivity unchanged."""
    if element_type in PROTEUS_HEX27_NAMES:
        return connectivity[:, proteus_hex27_to_exodus_hex27], "HEX27"
    return connectivity, element_type
