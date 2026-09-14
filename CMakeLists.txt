################################################################################
#
#  OBB_alg - OBB (Oriented Bounding Box) Tree algorithm built on top of CHAI3D.
#
#  Layout:
#    chai3d-master/chai3d-master/   original CHAI3D library (left unmodified)
#    Algorithm/                     OBB algorithm code (target: obb)
#
################################################################################

cmake_minimum_required (VERSION 3.19)
project (OBB_alg C CXX)

# Location of the CHAI3D source tree.
set (CHAI3D_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/chai3d-master/chai3d-master)

# CHAI3D declares cmake_minimum_required(VERSION 3.0.0), which CMake 4.x no
# longer accepts. Raise the policy floor for the CHAI3D subtree only.
set (CMAKE_POLICY_VERSION_MINIMUM 3.5)

# CHAI3D static library (target: chai3d). EXCLUDE_FROM_ALL keeps its examples
# and utilities out of the default build; chai3d itself is still built when a
# target links to it.
add_subdirectory (${CHAI3D_ROOT} chai3d EXCLUDE_FROM_ALL)

unset (CMAKE_POLICY_VERSION_MINIMUM)

# OBB algorithm library and its tests.
enable_testing ()
add_subdirectory (Algorithm)
