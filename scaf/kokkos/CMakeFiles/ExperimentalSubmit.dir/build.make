# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.28

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.28.1/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.28.1/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/hyw/Desktop/ERM/scaf

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/hyw/Desktop/ERM/scaf/scaf

# Utility rule file for ExperimentalSubmit.

# Include any custom commands dependencies for this target.
include kokkos/CMakeFiles/ExperimentalSubmit.dir/compiler_depend.make

# Include the progress variables for this target.
include kokkos/CMakeFiles/ExperimentalSubmit.dir/progress.make

kokkos/CMakeFiles/ExperimentalSubmit:
	cd /Users/hyw/Desktop/ERM/scaf/scaf/kokkos && /opt/homebrew/Cellar/cmake/3.28.1/bin/ctest -D ExperimentalSubmit

ExperimentalSubmit: kokkos/CMakeFiles/ExperimentalSubmit
ExperimentalSubmit: kokkos/CMakeFiles/ExperimentalSubmit.dir/build.make
.PHONY : ExperimentalSubmit

# Rule to build all files generated by this target.
kokkos/CMakeFiles/ExperimentalSubmit.dir/build: ExperimentalSubmit
.PHONY : kokkos/CMakeFiles/ExperimentalSubmit.dir/build

kokkos/CMakeFiles/ExperimentalSubmit.dir/clean:
	cd /Users/hyw/Desktop/ERM/scaf/scaf/kokkos && $(CMAKE_COMMAND) -P CMakeFiles/ExperimentalSubmit.dir/cmake_clean.cmake
.PHONY : kokkos/CMakeFiles/ExperimentalSubmit.dir/clean

kokkos/CMakeFiles/ExperimentalSubmit.dir/depend:
	cd /Users/hyw/Desktop/ERM/scaf/scaf && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/hyw/Desktop/ERM/scaf /Users/hyw/Desktop/ERM/scaf/kokkos /Users/hyw/Desktop/ERM/scaf/scaf /Users/hyw/Desktop/ERM/scaf/scaf/kokkos /Users/hyw/Desktop/ERM/scaf/scaf/kokkos/CMakeFiles/ExperimentalSubmit.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : kokkos/CMakeFiles/ExperimentalSubmit.dir/depend

