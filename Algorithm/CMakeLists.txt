################################################################################
#
#  OBB algorithm library (target: obb).
#
#  Linking to "obb" gives a target everything needed to use both the OBB code
#  and CHAI3D: include paths, compile definitions and libraries.
#
################################################################################

file (GLOB_RECURSE OBB_HEADERS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h)
file (GLOB_RECURSE OBB_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

# Header-only while there are no .cpp files yet; becomes a static library as
# soon as the first .cpp is added to src/.
if (OBB_SOURCES)
    add_library (obb STATIC ${OBB_HEADERS} ${OBB_SOURCES})
    set (OBB_SCOPE PUBLIC)
else ()
    add_library (obb INTERFACE ${OBB_HEADERS})
    set (OBB_SCOPE INTERFACE)
endif ()

# CHAI3D only sets its include paths and definitions at directory level, so
# they are not inherited through target_link_libraries and are repeated here.
target_include_directories (obb ${OBB_SCOPE}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CHAI3D_ROOT}/src
    ${CHAI3D_ROOT}/externals/Eigen
    ${CHAI3D_ROOT}/externals/glew/include)

target_compile_definitions (obb ${OBB_SCOPE} GLEW_STATIC)

find_package (OpenGL REQUIRED)
target_link_libraries (obb ${OBB_SCOPE} chai3d OpenGL::GL OpenGL::GLU)

if (WIN32)
    if (CMAKE_SIZEOF_VOID_P EQUAL 8)
        target_compile_definitions (obb ${OBB_SCOPE} WIN64)
    else ()
        target_compile_definitions (obb ${OBB_SCOPE} WIN32)
    endif ()
    if (MSVC)
        target_compile_definitions (obb ${OBB_SCOPE} _CRT_SECURE_NO_DEPRECATE)
    endif ()
    target_link_libraries (obb ${OBB_SCOPE} winmm)
endif ()

# Viewer: displays the OBB trees of a hand and an object (GLFW window, from CHAI3D extras)
if (TARGET glfw)
    add_executable (obb_viewer ${CMAKE_CURRENT_SOURCE_DIR}/examples/obb_viewer/main.cpp)
    target_include_directories (obb_viewer PRIVATE ${CHAI3D_ROOT}/extras/GLFW/include)
    target_link_libraries (obb_viewer PRIVATE obb glfw)
endif ()

# Benchmark: OBB tree vs CHAI3D AABB and exact references, 1k to 1M triangles
add_executable (obb_benchmark ${CMAKE_CURRENT_SOURCE_DIR}/examples/obb_benchmark/main.cpp)
target_link_libraries (obb_benchmark PRIVATE obb)

# Tests: one executable per tests/*.cpp file
file (GLOB OBB_TESTS CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
foreach (test_source ${OBB_TESTS})
    get_filename_component (test_name ${test_source} NAME_WE)
    add_executable (${test_name} ${test_source} ${CMAKE_CURRENT_SOURCE_DIR}/tests/TestUtils.h)
    target_link_libraries (${test_name} PRIVATE obb)
    add_test (NAME ${test_name} COMMAND ${test_name})
endforeach ()
