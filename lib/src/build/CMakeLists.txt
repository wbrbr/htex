include(CMakePackageConfigHelpers)

set(CMAKE_DIR "share/cmake/Htex")

write_basic_package_version_file("htex-config-version.cmake"
    VERSION "${HTEX_MAJOR_VERSION}.${HTEX_MINOR_VERSION}"
    COMPATIBILITY AnyNewerVersion)

configure_package_config_file(
    "htex-config.cmake" "htex-config.cmake"
    INSTALL_DESTINATION "${CMAKE_DIR}" PATH_VARS
    CMAKE_INSTALL_PREFIX CMAKE_INSTALL_LIBDIR CMAKE_INSTALL_INCLUDEDIR)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/htex-config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/htex-config-version.cmake"
    DESTINATION "${CMAKE_DIR}" COMPONENT "devel")

install(EXPORT Htex NAMESPACE Htex::
    FILE "htex-exports.cmake" DESTINATION "${CMAKE_DIR}" COMPONENT "devel")
