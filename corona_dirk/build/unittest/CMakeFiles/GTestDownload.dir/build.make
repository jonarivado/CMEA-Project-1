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
CMAKE_SOURCE_DIR = /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build

# Utility rule file for GTestDownload.

# Include any custom commands dependencies for this target.
include unittest/CMakeFiles/GTestDownload.dir/compiler_depend.make

# Include the progress variables for this target.
include unittest/CMakeFiles/GTestDownload.dir/progress.make

unittest/CMakeFiles/GTestDownload: unittest/CMakeFiles/GTestDownload-complete

unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-install
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-mkdir
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-download
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-update
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-patch
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-configure
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-build
unittest/CMakeFiles/GTestDownload-complete: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-install
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Completed 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/CMakeFiles
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/CMakeFiles/GTestDownload-complete
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-done

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-build: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-configure
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Performing build step for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary && $(MAKE)
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-build

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-configure: unittest/GTestDownload-prefix/tmp/GTestDownload-cfgcmd.txt
unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-configure: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-patch
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Performing configure step for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary && /usr/bin/cmake.exe -Dgtest_force_shared_crt=ON -DCMAKE_CXX_COMPILER=/usr/bin/c++.exe -DCMAKE_C_COMPILER=/usr/bin/cc -DCMAKE_INSTALL_PREFIX=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary "-GUnix Makefiles" /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_source
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-configure

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-download: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-urlinfo.txt
unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-download: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-mkdir
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Performing download step (download, verify and extract) for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -P /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/download-GTestDownload.cmake
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -P /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/verify-GTestDownload.cmake
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -P /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/extract-GTestDownload.cmake
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-download

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-install: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-build
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "No install step for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary && /usr/bin/cmake.exe -E echo_append
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-install

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-mkdir:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Creating directories for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_source
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/gtest_binary
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/tmp
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E make_directory /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-mkdir

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-patch: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-update
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "No patch step for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E echo_append
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-patch

unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-update: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-download
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "No update step for 'GTestDownload'"
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E echo_append
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && /usr/bin/cmake.exe -E touch /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-update

GTestDownload: unittest/CMakeFiles/GTestDownload
GTestDownload: unittest/CMakeFiles/GTestDownload-complete
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-build
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-configure
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-download
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-install
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-mkdir
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-patch
GTestDownload: unittest/GTestDownload-prefix/src/GTestDownload-stamp/GTestDownload-update
GTestDownload: unittest/CMakeFiles/GTestDownload.dir/build.make
.PHONY : GTestDownload

# Rule to build all files generated by this target.
unittest/CMakeFiles/GTestDownload.dir/build: GTestDownload
.PHONY : unittest/CMakeFiles/GTestDownload.dir/build

unittest/CMakeFiles/GTestDownload.dir/clean:
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest && $(CMAKE_COMMAND) -P CMakeFiles/GTestDownload.dir/cmake_clean.cmake
.PHONY : unittest/CMakeFiles/GTestDownload.dir/clean

unittest/CMakeFiles/GTestDownload.dir/depend:
	cd /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/unittest /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest /cygdrive/c/Users/jonat/OneDrive/Documents/CMEA/project1/corona_dirk/build/unittest/CMakeFiles/GTestDownload.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : unittest/CMakeFiles/GTestDownload.dir/depend

