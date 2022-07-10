configure_file(HtexVersion.h.in
    ${CMAKE_CURRENT_SOURCE_DIR}/HtexVersion.h @ONLY)

set(SRCS
    HtexCache.cpp
    HtexFilters.cpp
    HtexHalf.cpp
    HtexReader.cpp
    HtexSeparableFilter.cpp
    HtexSeparableKernel.cpp
    HtexTriangleFilter.cpp
    HtexTriangleKernel.cpp
    HtexUtils.cpp
    HtexWriter.cpp)

if(HTEX_BUILD_STATIC_LIBS)
    add_library(Htex_static STATIC ${SRCS})
    set_target_properties(Htex_static PROPERTIES OUTPUT_NAME Ptex)
    target_compile_definitions(Htex_static PUBLIC HTEX_STATIC)
    target_include_directories(Htex_static
    PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(Htex_static
        PUBLIC Threads::Threads PkgConfig::Ptex_ZLIB)
    install(TARGETS Htex_static EXPORT Htex DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(HTEX_BUILD_SHARED_LIBS)
    add_library(Htex_dynamic SHARED ${SRCS})
    set_target_properties(Htex_dynamic PROPERTIES
        OUTPUT_NAME Htex
        SOVERSION "${HTEX_MAJOR_VERSION}.${HTEX_MINOR_VERSION}")
    target_include_directories(Htex_dynamic
        PUBLIC
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(Htex_dynamic PRIVATE HTEX_EXPORTS)
    target_link_libraries(Htex_dynamic
        PUBLIC Threads::Threads PkgConfig::Ptex_ZLIB)
    install(TARGETS Htex_dynamic EXPORT Htex DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

install(FILES
        PtexExports.h
        PtexHalf.h
        PtexInt.h
        PtexPlatform.h
        Ptexture.h
        HtexUtils.h
        PtexVersion.h
        CatmullClark.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        COMPONENT devel)
