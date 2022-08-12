cmake_minimum_required (VERSION 3.1)

project(htex)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_AO_BAKER "Build the AO baker app (requires Embree)")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

if(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
  # Update if necessary
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "-g")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")

# disable GLFW docs, examples and tests
# see http://www.glfw.org/docs/latest/build_guide.html
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

set(PTEX_BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
set(PTEX_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(submodules/ptex)
set_property(TARGET Ptex_static ptxinfo halftest ftest rtest wtest PROPERTY FOLDER "ptex")
include_directories(submodules/ptex/src/ptex)

set(HTEX_BUILD_STATIC_LIBS ON CACHE BOOL "" FORCE)
set(HTEX_BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(submodules/htex-api)
set_property(TARGET Htex_static PROPERTY FOLDER "htex")
include_directories(submodules/htex-api/src/htex)


# set path to dependencies
add_subdirectory(submodules/glfw glfw)
include_directories(submodules/glfw/include)
include_directories(submodules/imgui)
include_directories(submodules/dj_opengl)
include_directories(submodules/stb)
include_directories(submodules/glm)
include_directories(submodules/HalfedgeCatmullClark)

if (BUILD_AO_BAKER)
    set(EMBREE_ISPC_SUPPORT OFF CACHE BOOL "" FORCE)
    add_subdirectory(submodules/embree)
    include_directories(submodules/embree/include)
endif()

# imgui source files
set(IMGUI_SRC_DIR submodules/imgui)
aux_source_directory(${IMGUI_SRC_DIR} IMGUI_SRC_FILES)

# Add dependencies
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# set path to common code
include_directories(common)

unset(SRC_FILES)
set(DEMO ptx2htx)
set(SRC_DIR apps/ptx2htx)
aux_source_directory(${SRC_DIR} SRC_FILES)
add_executable(${DEMO} ${SRC_FILES})
target_include_directories(${DEMO} PRIVATE ${SRC_DIR})
target_link_libraries(${DEMO} glfw Ptex_static Htex_static)
target_compile_definitions(
        ${DEMO} PUBLIC
        -DPATH_TO_SRC_DIRECTORY="${CMAKE_SOURCE_DIR}/${SRC_DIR}/"
)

unset(SRC_FILES)
set(DEMO htxview)
set(SRC_DIR apps/htxview)
aux_source_directory(${SRC_DIR} SRC_FILES)
add_executable(${DEMO} ${IMGUI_SRC_FILES} ${SRC_FILES} common/glad/glad.c common/imgui_impl.cpp)
target_include_directories(${DEMO} PRIVATE ${SRC_DIR})
target_link_libraries(${DEMO} glfw Htex_static)
target_compile_definitions(
        ${DEMO} PUBLIC
        -DPATH_TO_SRC_DIRECTORY="${CMAKE_SOURCE_DIR}/${SRC_DIR}/"
        -DPATH_TO_ROOT_DIRECTORY="${CMAKE_SOURCE_DIR}/"
        -DIMGUI_IMPL_OPENGL_LOADER_GLAD
)

unset(SRC_FILES)
set(DEMO painter)
set(SRC_DIR apps/painter)
aux_source_directory(${SRC_DIR} SRC_FILES)
add_executable(${DEMO} ${IMGUI_SRC_FILES} ${SRC_FILES} common/glad/glad.c common/imgui_impl.cpp)
target_include_directories(${DEMO} PRIVATE ${SRC_DIR})
target_link_libraries(${DEMO} glfw Htex_static)
target_compile_definitions(
        ${DEMO} PUBLIC
        -DPATH_TO_SRC_DIRECTORY="${CMAKE_SOURCE_DIR}/${SRC_DIR}/"
        -DPATH_TO_ROOT_DIRECTORY="${CMAKE_SOURCE_DIR}/"
        -DIMGUI_IMPL_OPENGL_LOADER_GLAD
)

unset(SRC_FILES)
set(DEMO uv2htx)
set(SRC_DIR apps/uv2htx)
aux_source_directory(${SRC_DIR} SRC_FILES)
add_executable(${DEMO} ${SRC_FILES})
target_include_directories(${DEMO} PRIVATE ${SRC_DIR})
target_link_libraries(${DEMO} Htex_static)

if(BUILD_AO_BAKER)
    unset(SRC_FILES)
    set(DEMO aobaker)
    set(SRC_DIR apps/aobaker)
    aux_source_directory(${SRC_DIR} SRC_FILES)
    add_executable(${DEMO} ${SRC_FILES})
    target_include_directories(${DEMO} PRIVATE ${SRC_DIR})
    target_link_libraries(${DEMO} Htex_static embree)
endif()

unset(DEMO)
unset(SRC_FILES)

add_executable(obj_to_ccm submodules/HalfedgeCatmullClark/examples/obj_to_ccm.c)