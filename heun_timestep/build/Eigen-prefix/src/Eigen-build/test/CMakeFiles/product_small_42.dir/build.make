# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

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
CMAKE_COMMAND = /usr/bin/cmake.exe

# The command to remove a file.
RM = /usr/bin/cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build

# Include any dependencies generated for this target.
include test/CMakeFiles/product_small_42.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include test/CMakeFiles/product_small_42.dir/compiler_depend.make

# Include the progress variables for this target.
include test/CMakeFiles/product_small_42.dir/progress.make

# Include the compile flags for this target's objects.
include test/CMakeFiles/product_small_42.dir/flags.make

test/CMakeFiles/product_small_42.dir/product_small.cpp.o: test/CMakeFiles/product_small_42.dir/flags.make
test/CMakeFiles/product_small_42.dir/product_small.cpp.o: /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/test/product_small.cpp
test/CMakeFiles/product_small_42.dir/product_small.cpp.o: test/CMakeFiles/product_small_42.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object test/CMakeFiles/product_small_42.dir/product_small.cpp.o"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test && /usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT test/CMakeFiles/product_small_42.dir/product_small.cpp.o -MF CMakeFiles/product_small_42.dir/product_small.cpp.o.d -o CMakeFiles/product_small_42.dir/product_small.cpp.o -c /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/test/product_small.cpp

test/CMakeFiles/product_small_42.dir/product_small.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/product_small_42.dir/product_small.cpp.i"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test && /usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/test/product_small.cpp > CMakeFiles/product_small_42.dir/product_small.cpp.i

test/CMakeFiles/product_small_42.dir/product_small.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/product_small_42.dir/product_small.cpp.s"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test && /usr/bin/c++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/test/product_small.cpp -o CMakeFiles/product_small_42.dir/product_small.cpp.s

# Object files for target product_small_42
product_small_42_OBJECTS = \
"CMakeFiles/product_small_42.dir/product_small.cpp.o"

# External object files for target product_small_42
product_small_42_EXTERNAL_OBJECTS =

test/product_small_42.exe: test/CMakeFiles/product_small_42.dir/product_small.cpp.o
test/product_small_42.exe: test/CMakeFiles/product_small_42.dir/build.make
test/product_small_42.exe: test/CMakeFiles/product_small_42.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable product_small_42.exe"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/product_small_42.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
test/CMakeFiles/product_small_42.dir/build: test/product_small_42.exe
.PHONY : test/CMakeFiles/product_small_42.dir/build

test/CMakeFiles/product_small_42.dir/clean:
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test && $(CMAKE_COMMAND) -P CMakeFiles/product_small_42.dir/cmake_clean.cmake
.PHONY : test/CMakeFiles/product_small_42.dir/clean

test/CMakeFiles/product_small_42.dir/depend:
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen/test /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/heun_timestep/build/Eigen-prefix/src/Eigen-build/test/CMakeFiles/product_small_42.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : test/CMakeFiles/product_small_42.dir/depend

