# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.14

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /data/desertfinch/speidel/relate

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /data/desertfinch/speidel/relate/build

# Include any dependencies generated for this target.
include include/src/CMakeFiles/relateShared.dir/depend.make

# Include the progress variables for this target.
include include/src/CMakeFiles/relateShared.dir/progress.make

# Include the compile flags for this target's objects.
include include/src/CMakeFiles/relateShared.dir/flags.make

include/src/CMakeFiles/relateShared.dir/filesystem.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/filesystem.cpp.o: ../include/src/filesystem.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object include/src/CMakeFiles/relateShared.dir/filesystem.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/filesystem.cpp.o -c /data/desertfinch/speidel/relate/include/src/filesystem.cpp

include/src/CMakeFiles/relateShared.dir/filesystem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/filesystem.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/filesystem.cpp > CMakeFiles/relateShared.dir/filesystem.cpp.i

include/src/CMakeFiles/relateShared.dir/filesystem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/filesystem.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/filesystem.cpp -o CMakeFiles/relateShared.dir/filesystem.cpp.s

include/src/CMakeFiles/relateShared.dir/plot.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/plot.cpp.o: ../include/src/plot.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object include/src/CMakeFiles/relateShared.dir/plot.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/plot.cpp.o -c /data/desertfinch/speidel/relate/include/src/plot.cpp

include/src/CMakeFiles/relateShared.dir/plot.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/plot.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/plot.cpp > CMakeFiles/relateShared.dir/plot.cpp.i

include/src/CMakeFiles/relateShared.dir/plot.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/plot.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/plot.cpp -o CMakeFiles/relateShared.dir/plot.cpp.s

include/src/CMakeFiles/relateShared.dir/fast_log.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/fast_log.cpp.o: ../include/src/fast_log.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object include/src/CMakeFiles/relateShared.dir/fast_log.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/fast_log.cpp.o -c /data/desertfinch/speidel/relate/include/src/fast_log.cpp

include/src/CMakeFiles/relateShared.dir/fast_log.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/fast_log.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/fast_log.cpp > CMakeFiles/relateShared.dir/fast_log.cpp.i

include/src/CMakeFiles/relateShared.dir/fast_log.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/fast_log.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/fast_log.cpp -o CMakeFiles/relateShared.dir/fast_log.cpp.s

include/src/CMakeFiles/relateShared.dir/collapsed_matrix.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/collapsed_matrix.cpp.o: ../include/src/collapsed_matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object include/src/CMakeFiles/relateShared.dir/collapsed_matrix.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/collapsed_matrix.cpp.o -c /data/desertfinch/speidel/relate/include/src/collapsed_matrix.cpp

include/src/CMakeFiles/relateShared.dir/collapsed_matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/collapsed_matrix.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/collapsed_matrix.cpp > CMakeFiles/relateShared.dir/collapsed_matrix.cpp.i

include/src/CMakeFiles/relateShared.dir/collapsed_matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/collapsed_matrix.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/collapsed_matrix.cpp -o CMakeFiles/relateShared.dir/collapsed_matrix.cpp.s

include/src/CMakeFiles/relateShared.dir/data.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/data.cpp.o: ../include/src/data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object include/src/CMakeFiles/relateShared.dir/data.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/data.cpp.o -c /data/desertfinch/speidel/relate/include/src/data.cpp

include/src/CMakeFiles/relateShared.dir/data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/data.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/data.cpp > CMakeFiles/relateShared.dir/data.cpp.i

include/src/CMakeFiles/relateShared.dir/data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/data.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/data.cpp -o CMakeFiles/relateShared.dir/data.cpp.s

include/src/CMakeFiles/relateShared.dir/sample.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/sample.cpp.o: ../include/src/sample.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object include/src/CMakeFiles/relateShared.dir/sample.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/sample.cpp.o -c /data/desertfinch/speidel/relate/include/src/sample.cpp

include/src/CMakeFiles/relateShared.dir/sample.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/sample.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/sample.cpp > CMakeFiles/relateShared.dir/sample.cpp.i

include/src/CMakeFiles/relateShared.dir/sample.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/sample.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/sample.cpp -o CMakeFiles/relateShared.dir/sample.cpp.s

include/src/CMakeFiles/relateShared.dir/fast_painting.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/fast_painting.cpp.o: ../include/src/fast_painting.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object include/src/CMakeFiles/relateShared.dir/fast_painting.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/fast_painting.cpp.o -c /data/desertfinch/speidel/relate/include/src/fast_painting.cpp

include/src/CMakeFiles/relateShared.dir/fast_painting.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/fast_painting.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/fast_painting.cpp > CMakeFiles/relateShared.dir/fast_painting.cpp.i

include/src/CMakeFiles/relateShared.dir/fast_painting.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/fast_painting.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/fast_painting.cpp -o CMakeFiles/relateShared.dir/fast_painting.cpp.s

include/src/CMakeFiles/relateShared.dir/anc.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/anc.cpp.o: ../include/src/anc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object include/src/CMakeFiles/relateShared.dir/anc.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/anc.cpp.o -c /data/desertfinch/speidel/relate/include/src/anc.cpp

include/src/CMakeFiles/relateShared.dir/anc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/anc.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/anc.cpp > CMakeFiles/relateShared.dir/anc.cpp.i

include/src/CMakeFiles/relateShared.dir/anc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/anc.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/anc.cpp -o CMakeFiles/relateShared.dir/anc.cpp.s

include/src/CMakeFiles/relateShared.dir/mutations.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/mutations.cpp.o: ../include/src/mutations.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object include/src/CMakeFiles/relateShared.dir/mutations.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/mutations.cpp.o -c /data/desertfinch/speidel/relate/include/src/mutations.cpp

include/src/CMakeFiles/relateShared.dir/mutations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/mutations.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/mutations.cpp > CMakeFiles/relateShared.dir/mutations.cpp.i

include/src/CMakeFiles/relateShared.dir/mutations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/mutations.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/mutations.cpp -o CMakeFiles/relateShared.dir/mutations.cpp.s

include/src/CMakeFiles/relateShared.dir/tree_builder.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/tree_builder.cpp.o: ../include/src/tree_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object include/src/CMakeFiles/relateShared.dir/tree_builder.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/tree_builder.cpp.o -c /data/desertfinch/speidel/relate/include/src/tree_builder.cpp

include/src/CMakeFiles/relateShared.dir/tree_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/tree_builder.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/tree_builder.cpp > CMakeFiles/relateShared.dir/tree_builder.cpp.i

include/src/CMakeFiles/relateShared.dir/tree_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/tree_builder.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/tree_builder.cpp -o CMakeFiles/relateShared.dir/tree_builder.cpp.s

include/src/CMakeFiles/relateShared.dir/anc_builder.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/anc_builder.cpp.o: ../include/src/anc_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object include/src/CMakeFiles/relateShared.dir/anc_builder.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/anc_builder.cpp.o -c /data/desertfinch/speidel/relate/include/src/anc_builder.cpp

include/src/CMakeFiles/relateShared.dir/anc_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/anc_builder.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/anc_builder.cpp > CMakeFiles/relateShared.dir/anc_builder.cpp.i

include/src/CMakeFiles/relateShared.dir/anc_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/anc_builder.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/anc_builder.cpp -o CMakeFiles/relateShared.dir/anc_builder.cpp.s

include/src/CMakeFiles/relateShared.dir/tree_comparer.cpp.o: include/src/CMakeFiles/relateShared.dir/flags.make
include/src/CMakeFiles/relateShared.dir/tree_comparer.cpp.o: ../include/src/tree_comparer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object include/src/CMakeFiles/relateShared.dir/tree_comparer.cpp.o"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateShared.dir/tree_comparer.cpp.o -c /data/desertfinch/speidel/relate/include/src/tree_comparer.cpp

include/src/CMakeFiles/relateShared.dir/tree_comparer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateShared.dir/tree_comparer.cpp.i"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /data/desertfinch/speidel/relate/include/src/tree_comparer.cpp > CMakeFiles/relateShared.dir/tree_comparer.cpp.i

include/src/CMakeFiles/relateShared.dir/tree_comparer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateShared.dir/tree_comparer.cpp.s"
	cd /data/desertfinch/speidel/relate/build/include/src && /bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /data/desertfinch/speidel/relate/include/src/tree_comparer.cpp -o CMakeFiles/relateShared.dir/tree_comparer.cpp.s

# Object files for target relateShared
relateShared_OBJECTS = \
"CMakeFiles/relateShared.dir/filesystem.cpp.o" \
"CMakeFiles/relateShared.dir/plot.cpp.o" \
"CMakeFiles/relateShared.dir/fast_log.cpp.o" \
"CMakeFiles/relateShared.dir/collapsed_matrix.cpp.o" \
"CMakeFiles/relateShared.dir/data.cpp.o" \
"CMakeFiles/relateShared.dir/sample.cpp.o" \
"CMakeFiles/relateShared.dir/fast_painting.cpp.o" \
"CMakeFiles/relateShared.dir/anc.cpp.o" \
"CMakeFiles/relateShared.dir/mutations.cpp.o" \
"CMakeFiles/relateShared.dir/tree_builder.cpp.o" \
"CMakeFiles/relateShared.dir/anc_builder.cpp.o" \
"CMakeFiles/relateShared.dir/tree_comparer.cpp.o"

# External object files for target relateShared
relateShared_EXTERNAL_OBJECTS =

../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/filesystem.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/plot.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/fast_log.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/collapsed_matrix.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/data.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/sample.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/fast_painting.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/anc.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/mutations.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/tree_builder.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/anc_builder.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/tree_comparer.cpp.o
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/build.make
../bin/librelateShared.so: ../bin/libgzstreamShared.so
../bin/librelateShared.so: /usr/lib64/libz.so
../bin/librelateShared.so: include/src/CMakeFiles/relateShared.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/data/desertfinch/speidel/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX shared library ../../../bin/librelateShared.so"
	cd /data/desertfinch/speidel/relate/build/include/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/relateShared.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
include/src/CMakeFiles/relateShared.dir/build: ../bin/librelateShared.so

.PHONY : include/src/CMakeFiles/relateShared.dir/build

include/src/CMakeFiles/relateShared.dir/clean:
	cd /data/desertfinch/speidel/relate/build/include/src && $(CMAKE_COMMAND) -P CMakeFiles/relateShared.dir/cmake_clean.cmake
.PHONY : include/src/CMakeFiles/relateShared.dir/clean

include/src/CMakeFiles/relateShared.dir/depend:
	cd /data/desertfinch/speidel/relate/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /data/desertfinch/speidel/relate /data/desertfinch/speidel/relate/include/src /data/desertfinch/speidel/relate/build /data/desertfinch/speidel/relate/build/include/src /data/desertfinch/speidel/relate/build/include/src/CMakeFiles/relateShared.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : include/src/CMakeFiles/relateShared.dir/depend

