# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_COMMAND = /gpfs0/apps/well/cmake/3.7.2/bin/cmake

# The command to remove a file.
RM = /gpfs0/apps/well/cmake/3.7.2/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /users/myers/speidel/Documents/genomics/relate

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /users/myers/speidel/Documents/genomics/relate/build

# Include any dependencies generated for this target.
include include/src/CMakeFiles/relateStatic.dir/depend.make

# Include the progress variables for this target.
include include/src/CMakeFiles/relateStatic.dir/progress.make

# Include the compile flags for this target's objects.
include include/src/CMakeFiles/relateStatic.dir/flags.make

include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o: ../include/src/filesystem.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/filesystem.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/filesystem.cpp

include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/filesystem.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/filesystem.cpp > CMakeFiles/relateStatic.dir/filesystem.cpp.i

include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/filesystem.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/filesystem.cpp -o CMakeFiles/relateStatic.dir/filesystem.cpp.s

include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o


include/src/CMakeFiles/relateStatic.dir/plot.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/plot.cpp.o: ../include/src/plot.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object include/src/CMakeFiles/relateStatic.dir/plot.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/plot.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/plot.cpp

include/src/CMakeFiles/relateStatic.dir/plot.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/plot.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/plot.cpp > CMakeFiles/relateStatic.dir/plot.cpp.i

include/src/CMakeFiles/relateStatic.dir/plot.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/plot.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/plot.cpp -o CMakeFiles/relateStatic.dir/plot.cpp.s

include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/plot.cpp.o


include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o: ../include/src/fast_log.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/fast_log.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/fast_log.cpp

include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/fast_log.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/fast_log.cpp > CMakeFiles/relateStatic.dir/fast_log.cpp.i

include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/fast_log.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/fast_log.cpp -o CMakeFiles/relateStatic.dir/fast_log.cpp.s

include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o


include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o: ../include/src/collapsed_matrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/collapsed_matrix.cpp

include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/collapsed_matrix.cpp > CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.i

include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/collapsed_matrix.cpp -o CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.s

include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o


include/src/CMakeFiles/relateStatic.dir/data.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/data.cpp.o: ../include/src/data.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object include/src/CMakeFiles/relateStatic.dir/data.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/data.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/data.cpp

include/src/CMakeFiles/relateStatic.dir/data.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/data.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/data.cpp > CMakeFiles/relateStatic.dir/data.cpp.i

include/src/CMakeFiles/relateStatic.dir/data.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/data.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/data.cpp -o CMakeFiles/relateStatic.dir/data.cpp.s

include/src/CMakeFiles/relateStatic.dir/data.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/data.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/data.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/data.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/data.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/data.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/data.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/data.cpp.o


include/src/CMakeFiles/relateStatic.dir/sample.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/sample.cpp.o: ../include/src/sample.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object include/src/CMakeFiles/relateStatic.dir/sample.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/sample.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/sample.cpp

include/src/CMakeFiles/relateStatic.dir/sample.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/sample.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/sample.cpp > CMakeFiles/relateStatic.dir/sample.cpp.i

include/src/CMakeFiles/relateStatic.dir/sample.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/sample.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/sample.cpp -o CMakeFiles/relateStatic.dir/sample.cpp.s

include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/sample.cpp.o


include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o: ../include/src/fast_painting.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/fast_painting.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/fast_painting.cpp

include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/fast_painting.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/fast_painting.cpp > CMakeFiles/relateStatic.dir/fast_painting.cpp.i

include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/fast_painting.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/fast_painting.cpp -o CMakeFiles/relateStatic.dir/fast_painting.cpp.s

include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o


include/src/CMakeFiles/relateStatic.dir/anc.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/anc.cpp.o: ../include/src/anc.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object include/src/CMakeFiles/relateStatic.dir/anc.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/anc.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/anc.cpp

include/src/CMakeFiles/relateStatic.dir/anc.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/anc.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/anc.cpp > CMakeFiles/relateStatic.dir/anc.cpp.i

include/src/CMakeFiles/relateStatic.dir/anc.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/anc.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/anc.cpp -o CMakeFiles/relateStatic.dir/anc.cpp.s

include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/anc.cpp.o


include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o: ../include/src/mutations.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/mutations.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/mutations.cpp

include/src/CMakeFiles/relateStatic.dir/mutations.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/mutations.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/mutations.cpp > CMakeFiles/relateStatic.dir/mutations.cpp.i

include/src/CMakeFiles/relateStatic.dir/mutations.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/mutations.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/mutations.cpp -o CMakeFiles/relateStatic.dir/mutations.cpp.s

include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o


include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o: ../include/src/tree_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/tree_builder.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/tree_builder.cpp

include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/tree_builder.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/tree_builder.cpp > CMakeFiles/relateStatic.dir/tree_builder.cpp.i

include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/tree_builder.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/tree_builder.cpp -o CMakeFiles/relateStatic.dir/tree_builder.cpp.s

include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o


include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o: ../include/src/branch_length_estimator.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/branch_length_estimator.cpp

include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/branch_length_estimator.cpp > CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.i

include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/branch_length_estimator.cpp -o CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.s

include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o


include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o: ../include/src/anc_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/anc_builder.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/anc_builder.cpp

include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/anc_builder.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/anc_builder.cpp > CMakeFiles/relateStatic.dir/anc_builder.cpp.i

include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/anc_builder.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/anc_builder.cpp -o CMakeFiles/relateStatic.dir/anc_builder.cpp.s

include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o


include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o: include/src/CMakeFiles/relateStatic.dir/flags.make
include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o: ../include/src/tree_comparer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/relateStatic.dir/tree_comparer.cpp.o -c /users/myers/speidel/Documents/genomics/relate/include/src/tree_comparer.cpp

include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/relateStatic.dir/tree_comparer.cpp.i"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /users/myers/speidel/Documents/genomics/relate/include/src/tree_comparer.cpp > CMakeFiles/relateStatic.dir/tree_comparer.cpp.i

include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/relateStatic.dir/tree_comparer.cpp.s"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && /apps/well/gcc/5.4.0/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /users/myers/speidel/Documents/genomics/relate/include/src/tree_comparer.cpp -o CMakeFiles/relateStatic.dir/tree_comparer.cpp.s

include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.requires:

.PHONY : include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.requires

include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.provides: include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.requires
	$(MAKE) -f include/src/CMakeFiles/relateStatic.dir/build.make include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.provides.build
.PHONY : include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.provides

include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.provides.build: include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o


# Object files for target relateStatic
relateStatic_OBJECTS = \
"CMakeFiles/relateStatic.dir/filesystem.cpp.o" \
"CMakeFiles/relateStatic.dir/plot.cpp.o" \
"CMakeFiles/relateStatic.dir/fast_log.cpp.o" \
"CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o" \
"CMakeFiles/relateStatic.dir/data.cpp.o" \
"CMakeFiles/relateStatic.dir/sample.cpp.o" \
"CMakeFiles/relateStatic.dir/fast_painting.cpp.o" \
"CMakeFiles/relateStatic.dir/anc.cpp.o" \
"CMakeFiles/relateStatic.dir/mutations.cpp.o" \
"CMakeFiles/relateStatic.dir/tree_builder.cpp.o" \
"CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o" \
"CMakeFiles/relateStatic.dir/anc_builder.cpp.o" \
"CMakeFiles/relateStatic.dir/tree_comparer.cpp.o"

# External object files for target relateStatic
relateStatic_EXTERNAL_OBJECTS =

../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/plot.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/data.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/sample.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/anc.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/build.make
../bin/librelateStatic.a: include/src/CMakeFiles/relateStatic.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/users/myers/speidel/Documents/genomics/relate/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX static library ../../../bin/librelateStatic.a"
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && $(CMAKE_COMMAND) -P CMakeFiles/relateStatic.dir/cmake_clean_target.cmake
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/relateStatic.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
include/src/CMakeFiles/relateStatic.dir/build: ../bin/librelateStatic.a

.PHONY : include/src/CMakeFiles/relateStatic.dir/build

include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/filesystem.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/plot.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/fast_log.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/collapsed_matrix.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/data.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/sample.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/fast_painting.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/anc.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/mutations.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/tree_builder.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/branch_length_estimator.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/anc_builder.cpp.o.requires
include/src/CMakeFiles/relateStatic.dir/requires: include/src/CMakeFiles/relateStatic.dir/tree_comparer.cpp.o.requires

.PHONY : include/src/CMakeFiles/relateStatic.dir/requires

include/src/CMakeFiles/relateStatic.dir/clean:
	cd /users/myers/speidel/Documents/genomics/relate/build/include/src && $(CMAKE_COMMAND) -P CMakeFiles/relateStatic.dir/cmake_clean.cmake
.PHONY : include/src/CMakeFiles/relateStatic.dir/clean

include/src/CMakeFiles/relateStatic.dir/depend:
	cd /users/myers/speidel/Documents/genomics/relate/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /users/myers/speidel/Documents/genomics/relate /users/myers/speidel/Documents/genomics/relate/include/src /users/myers/speidel/Documents/genomics/relate/build /users/myers/speidel/Documents/genomics/relate/build/include/src /users/myers/speidel/Documents/genomics/relate/build/include/src/CMakeFiles/relateStatic.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : include/src/CMakeFiles/relateStatic.dir/depend

