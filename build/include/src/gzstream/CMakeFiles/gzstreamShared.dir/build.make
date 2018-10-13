# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/leo/Documents/genomics/relate

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/leo/Documents/genomics/relate/build

# Include any dependencies generated for this target.
include include/src/gzstream/CMakeFiles/gzstreamShared.dir/depend.make

# Include the progress variables for this target.
include include/src/gzstream/CMakeFiles/gzstreamShared.dir/progress.make

# Include the compile flags for this target's objects.
include include/src/gzstream/CMakeFiles/gzstreamShared.dir/flags.make

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o: include/src/gzstream/CMakeFiles/gzstreamShared.dir/flags.make
include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o: ../include/src/gzstream/gzstream.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/leo/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o"
	cd /Users/leo/Documents/genomics/relate/build/include/src/gzstream && /Library/Developer/CommandLineTools/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gzstreamShared.dir/gzstream.cpp.o -c /Users/leo/Documents/genomics/relate/include/src/gzstream/gzstream.cpp

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gzstreamShared.dir/gzstream.cpp.i"
	cd /Users/leo/Documents/genomics/relate/build/include/src/gzstream && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/leo/Documents/genomics/relate/include/src/gzstream/gzstream.cpp > CMakeFiles/gzstreamShared.dir/gzstream.cpp.i

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gzstreamShared.dir/gzstream.cpp.s"
	cd /Users/leo/Documents/genomics/relate/build/include/src/gzstream && /Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/leo/Documents/genomics/relate/include/src/gzstream/gzstream.cpp -o CMakeFiles/gzstreamShared.dir/gzstream.cpp.s

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.requires:

.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.requires

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.provides: include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.requires
	$(MAKE) -f include/src/gzstream/CMakeFiles/gzstreamShared.dir/build.make include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.provides.build
.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.provides

include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.provides.build: include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o


# Object files for target gzstreamShared
gzstreamShared_OBJECTS = \
"CMakeFiles/gzstreamShared.dir/gzstream.cpp.o"

# External object files for target gzstreamShared
gzstreamShared_EXTERNAL_OBJECTS =

../bin/libgzstreamShared.dylib: include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o
../bin/libgzstreamShared.dylib: include/src/gzstream/CMakeFiles/gzstreamShared.dir/build.make
../bin/libgzstreamShared.dylib: /opt/local/lib/libz.dylib
../bin/libgzstreamShared.dylib: include/src/gzstream/CMakeFiles/gzstreamShared.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/leo/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX shared library ../../../../bin/libgzstreamShared.dylib"
	cd /Users/leo/Documents/genomics/relate/build/include/src/gzstream && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gzstreamShared.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
include/src/gzstream/CMakeFiles/gzstreamShared.dir/build: ../bin/libgzstreamShared.dylib

.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/build

include/src/gzstream/CMakeFiles/gzstreamShared.dir/requires: include/src/gzstream/CMakeFiles/gzstreamShared.dir/gzstream.cpp.o.requires

.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/requires

include/src/gzstream/CMakeFiles/gzstreamShared.dir/clean:
	cd /Users/leo/Documents/genomics/relate/build/include/src/gzstream && $(CMAKE_COMMAND) -P CMakeFiles/gzstreamShared.dir/cmake_clean.cmake
.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/clean

include/src/gzstream/CMakeFiles/gzstreamShared.dir/depend:
	cd /Users/leo/Documents/genomics/relate/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/leo/Documents/genomics/relate /Users/leo/Documents/genomics/relate/include/src/gzstream /Users/leo/Documents/genomics/relate/build /Users/leo/Documents/genomics/relate/build/include/src/gzstream /Users/leo/Documents/genomics/relate/build/include/src/gzstream/CMakeFiles/gzstreamShared.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : include/src/gzstream/CMakeFiles/gzstreamShared.dir/depend

