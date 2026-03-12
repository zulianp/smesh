#!/usr/bin/env bash

set -e

# ./square TRI3 4 4 0 0 1 1 tri3_square
# mesh=tri3_square

mesh=pump
ranks=8

SMESH_ORDERING_TYPE=hilbert3 ./sfc "$mesh" "$mesh"
$CODE_DIR/smesh/python/smesh/raw_to_db.py "$mesh" "$mesh".vtk

mpiexec -np 2 ./skin "$mesh" "$mesh"_skin
$CODE_DIR/smesh/python/smesh/raw_to_db.py "$mesh"_skin "$mesh"_skin.vtk


rm -rf test_0 test_1

mpiexec -np $ranks ./io_test "$mesh" "$mesh"_dump

$CODE_DIR/smesh/python/smesh/raw_to_db.py test_0 test_0.vtk -c test_0/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_1 test_1.vtk -c test_1/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_2 test_2.vtk -c test_2/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_3 test_3.vtk -c test_3/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_4 test_4.vtk -c test_4/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_5 test_5.vtk -c test_5/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_6 test_6.vtk -c test_6/elemental_data.int64
$CODE_DIR/smesh/python/smesh/raw_to_db.py test_7 test_7.vtk -c test_7/elemental_data.int64