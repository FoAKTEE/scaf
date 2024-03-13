# Install script for directory: /Users/hyw/Desktop/ERM/scaf/kokkos

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/core/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/containers/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/algorithms/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/example/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos" TYPE FILE FILES
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosConfig.cmake"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosConfigCommon.cmake"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosConfigVersion.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets.cmake"
         "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/CMakeFiles/Export/7c5c9398a7e290ebb35584c8631204e4/KokkosTargets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos/KokkosTargets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos" TYPE FILE FILES "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/CMakeFiles/Export/7c5c9398a7e290ebb35584c8631204e4/KokkosTargets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Rr][Ee][Ll][Ee][Aa][Ss][Ee])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/Kokkos" TYPE FILE FILES "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/CMakeFiles/Export/7c5c9398a7e290ebb35584c8631204e4/KokkosTargets-release.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/kokkos" TYPE FILE FILES "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosCore_config.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE PROGRAM FILES
    "/Users/hyw/Desktop/ERM/scaf/kokkos/bin/nvcc_wrapper"
    "/Users/hyw/Desktop/ERM/scaf/kokkos/bin/hpcbind"
    "/Users/hyw/Desktop/ERM/scaf/kokkos/bin/kokkos_launch_compiler"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/temp/kokkos_launch_compiler"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/kokkos" TYPE FILE FILES
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosCore_config.h"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosCore_Config_FwdBackend.hpp"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosCore_Config_SetupBackend.hpp"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosCore_Config_DeclareBackend.hpp"
    "/Users/hyw/Desktop/ERM/scaf/scaf/kokkos/KokkosCore_Config_PostInclude.hpp"
    )
endif()
