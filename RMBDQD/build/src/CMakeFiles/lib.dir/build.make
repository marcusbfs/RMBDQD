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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/marcus/workspace/programa_com_malhas/RMBDQD

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/marcus/workspace/programa_com_malhas/RMBDQD/build

# Include any dependencies generated for this target.
include src/CMakeFiles/lib.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/lib.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/lib.dir/flags.make

src/CMakeFiles/lib.dir/mod/decomp.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/decomp.f90.o: ../src/mod/decomp.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building Fortran object src/CMakeFiles/lib.dir/mod/decomp.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/decomp.f90 -o CMakeFiles/lib.dir/mod/decomp.f90.o

src/CMakeFiles/lib.dir/mod/decomp.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/decomp.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/decomp.f90 > CMakeFiles/lib.dir/mod/decomp.f90.i

src/CMakeFiles/lib.dir/mod/decomp.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/decomp.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/decomp.f90 -o CMakeFiles/lib.dir/mod/decomp.f90.s

src/CMakeFiles/lib.dir/mod/decomp.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/decomp.f90.o.requires

src/CMakeFiles/lib.dir/mod/decomp.f90.o.provides: src/CMakeFiles/lib.dir/mod/decomp.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/decomp.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/decomp.f90.o.provides

src/CMakeFiles/lib.dir/mod/decomp.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/decomp.f90.o


src/CMakeFiles/lib.dir/mod/jacobi.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/jacobi.f90.o: ../src/mod/jacobi.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building Fortran object src/CMakeFiles/lib.dir/mod/jacobi.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/jacobi.f90 -o CMakeFiles/lib.dir/mod/jacobi.f90.o

src/CMakeFiles/lib.dir/mod/jacobi.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/jacobi.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/jacobi.f90 > CMakeFiles/lib.dir/mod/jacobi.f90.i

src/CMakeFiles/lib.dir/mod/jacobi.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/jacobi.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/jacobi.f90 -o CMakeFiles/lib.dir/mod/jacobi.f90.s

src/CMakeFiles/lib.dir/mod/jacobi.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/jacobi.f90.o.requires

src/CMakeFiles/lib.dir/mod/jacobi.f90.o.provides: src/CMakeFiles/lib.dir/mod/jacobi.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/jacobi.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/jacobi.f90.o.provides

src/CMakeFiles/lib.dir/mod/jacobi.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/jacobi.f90.o


src/CMakeFiles/lib.dir/mod/mult.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/mult.f90.o: ../src/mod/mult.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building Fortran object src/CMakeFiles/lib.dir/mod/mult.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/mult.f90 -o CMakeFiles/lib.dir/mod/mult.f90.o

src/CMakeFiles/lib.dir/mod/mult.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/mult.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/mult.f90 > CMakeFiles/lib.dir/mod/mult.f90.i

src/CMakeFiles/lib.dir/mod/mult.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/mult.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/mult.f90 -o CMakeFiles/lib.dir/mod/mult.f90.s

src/CMakeFiles/lib.dir/mod/mult.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/mult.f90.o.requires

src/CMakeFiles/lib.dir/mod/mult.f90.o.provides: src/CMakeFiles/lib.dir/mod/mult.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/mult.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/mult.f90.o.provides

src/CMakeFiles/lib.dir/mod/mult.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/mult.f90.o


src/CMakeFiles/lib.dir/mod/redbak.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/redbak.f90.o: ../src/mod/redbak.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building Fortran object src/CMakeFiles/lib.dir/mod/redbak.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/redbak.f90 -o CMakeFiles/lib.dir/mod/redbak.f90.o

src/CMakeFiles/lib.dir/mod/redbak.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/redbak.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/redbak.f90 > CMakeFiles/lib.dir/mod/redbak.f90.i

src/CMakeFiles/lib.dir/mod/redbak.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/redbak.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/redbak.f90 -o CMakeFiles/lib.dir/mod/redbak.f90.s

src/CMakeFiles/lib.dir/mod/redbak.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/redbak.f90.o.requires

src/CMakeFiles/lib.dir/mod/redbak.f90.o.provides: src/CMakeFiles/lib.dir/mod/redbak.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/redbak.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/redbak.f90.o.provides

src/CMakeFiles/lib.dir/mod/redbak.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/redbak.f90.o


src/CMakeFiles/lib.dir/mod/scheck.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/scheck.f90.o: ../src/mod/scheck.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building Fortran object src/CMakeFiles/lib.dir/mod/scheck.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/scheck.f90 -o CMakeFiles/lib.dir/mod/scheck.f90.o

src/CMakeFiles/lib.dir/mod/scheck.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/scheck.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/scheck.f90 > CMakeFiles/lib.dir/mod/scheck.f90.i

src/CMakeFiles/lib.dir/mod/scheck.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/scheck.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/scheck.f90 -o CMakeFiles/lib.dir/mod/scheck.f90.s

src/CMakeFiles/lib.dir/mod/scheck.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/scheck.f90.o.requires

src/CMakeFiles/lib.dir/mod/scheck.f90.o.provides: src/CMakeFiles/lib.dir/mod/scheck.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/scheck.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/scheck.f90.o.provides

src/CMakeFiles/lib.dir/mod/scheck.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/scheck.f90.o


src/CMakeFiles/lib.dir/mod/commonutils.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/commonutils.f90.o: ../src/mod/commonutils.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building Fortran object src/CMakeFiles/lib.dir/mod/commonutils.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/commonutils.f90 -o CMakeFiles/lib.dir/mod/commonutils.f90.o

src/CMakeFiles/lib.dir/mod/commonutils.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/commonutils.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/commonutils.f90 > CMakeFiles/lib.dir/mod/commonutils.f90.i

src/CMakeFiles/lib.dir/mod/commonutils.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/commonutils.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/commonutils.f90 -o CMakeFiles/lib.dir/mod/commonutils.f90.s

src/CMakeFiles/lib.dir/mod/commonutils.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/commonutils.f90.o.requires

src/CMakeFiles/lib.dir/mod/commonutils.f90.o.provides: src/CMakeFiles/lib.dir/mod/commonutils.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/commonutils.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/commonutils.f90.o.provides

src/CMakeFiles/lib.dir/mod/commonutils.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/commonutils.f90.o


src/CMakeFiles/lib.dir/mod/subrotinas.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/subrotinas.f90.o: ../src/mod/subrotinas.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building Fortran object src/CMakeFiles/lib.dir/mod/subrotinas.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/subrotinas.f90 -o CMakeFiles/lib.dir/mod/subrotinas.f90.o

src/CMakeFiles/lib.dir/mod/subrotinas.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/subrotinas.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/subrotinas.f90 > CMakeFiles/lib.dir/mod/subrotinas.f90.i

src/CMakeFiles/lib.dir/mod/subrotinas.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/subrotinas.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/subrotinas.f90 -o CMakeFiles/lib.dir/mod/subrotinas.f90.s

src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.requires

src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.provides: src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.provides

src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/subrotinas.f90.o


src/CMakeFiles/lib.dir/mod/subspace.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/subspace.f90.o: ../src/mod/subspace.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building Fortran object src/CMakeFiles/lib.dir/mod/subspace.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/subspace.f90 -o CMakeFiles/lib.dir/mod/subspace.f90.o

src/CMakeFiles/lib.dir/mod/subspace.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/subspace.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/subspace.f90 > CMakeFiles/lib.dir/mod/subspace.f90.i

src/CMakeFiles/lib.dir/mod/subspace.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/subspace.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/subspace.f90 -o CMakeFiles/lib.dir/mod/subspace.f90.s

src/CMakeFiles/lib.dir/mod/subspace.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/subspace.f90.o.requires

src/CMakeFiles/lib.dir/mod/subspace.f90.o.provides: src/CMakeFiles/lib.dir/mod/subspace.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/subspace.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/subspace.f90.o.provides

src/CMakeFiles/lib.dir/mod/subspace.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/subspace.f90.o


src/CMakeFiles/lib.dir/mod/help.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/help.f90.o: ../src/mod/help.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building Fortran object src/CMakeFiles/lib.dir/mod/help.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/help.f90 -o CMakeFiles/lib.dir/mod/help.f90.o

src/CMakeFiles/lib.dir/mod/help.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/help.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/help.f90 > CMakeFiles/lib.dir/mod/help.f90.i

src/CMakeFiles/lib.dir/mod/help.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/help.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/help.f90 -o CMakeFiles/lib.dir/mod/help.f90.s

src/CMakeFiles/lib.dir/mod/help.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/help.f90.o.requires

src/CMakeFiles/lib.dir/mod/help.f90.o.provides: src/CMakeFiles/lib.dir/mod/help.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/help.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/help.f90.o.provides

src/CMakeFiles/lib.dir/mod/help.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/help.f90.o


src/CMakeFiles/lib.dir/mod/f90getopt.f90.o: src/CMakeFiles/lib.dir/flags.make
src/CMakeFiles/lib.dir/mod/f90getopt.f90.o: ../src/mod/f90getopt.f90
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building Fortran object src/CMakeFiles/lib.dir/mod/f90getopt.f90.o"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -c /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/f90getopt.f90 -o CMakeFiles/lib.dir/mod/f90getopt.f90.o

src/CMakeFiles/lib.dir/mod/f90getopt.f90.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing Fortran source to CMakeFiles/lib.dir/mod/f90getopt.f90.i"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -E /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/f90getopt.f90 > CMakeFiles/lib.dir/mod/f90getopt.f90.i

src/CMakeFiles/lib.dir/mod/f90getopt.f90.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling Fortran source to assembly CMakeFiles/lib.dir/mod/f90getopt.f90.s"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && /usr/bin/gfortran  $(Fortran_DEFINES) $(Fortran_INCLUDES) $(Fortran_FLAGS) -S /home/marcus/workspace/programa_com_malhas/RMBDQD/src/mod/f90getopt.f90 -o CMakeFiles/lib.dir/mod/f90getopt.f90.s

src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.requires:

.PHONY : src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.requires

src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.provides: src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.requires
	$(MAKE) -f src/CMakeFiles/lib.dir/build.make src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.provides.build
.PHONY : src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.provides

src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.provides.build: src/CMakeFiles/lib.dir/mod/f90getopt.f90.o


# Object files for target lib
lib_OBJECTS = \
"CMakeFiles/lib.dir/mod/decomp.f90.o" \
"CMakeFiles/lib.dir/mod/jacobi.f90.o" \
"CMakeFiles/lib.dir/mod/mult.f90.o" \
"CMakeFiles/lib.dir/mod/redbak.f90.o" \
"CMakeFiles/lib.dir/mod/scheck.f90.o" \
"CMakeFiles/lib.dir/mod/commonutils.f90.o" \
"CMakeFiles/lib.dir/mod/subrotinas.f90.o" \
"CMakeFiles/lib.dir/mod/subspace.f90.o" \
"CMakeFiles/lib.dir/mod/help.f90.o" \
"CMakeFiles/lib.dir/mod/f90getopt.f90.o"

# External object files for target lib
lib_EXTERNAL_OBJECTS =

src/liblib.a: src/CMakeFiles/lib.dir/mod/decomp.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/jacobi.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/mult.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/redbak.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/scheck.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/commonutils.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/subrotinas.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/subspace.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/help.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/mod/f90getopt.f90.o
src/liblib.a: src/CMakeFiles/lib.dir/build.make
src/liblib.a: src/CMakeFiles/lib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marcus/workspace/programa_com_malhas/RMBDQD/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Linking Fortran static library liblib.a"
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && $(CMAKE_COMMAND) -P CMakeFiles/lib.dir/cmake_clean_target.cmake
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lib.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/lib.dir/build: src/liblib.a

.PHONY : src/CMakeFiles/lib.dir/build

src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/decomp.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/jacobi.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/mult.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/redbak.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/scheck.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/commonutils.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/subrotinas.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/subspace.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/help.f90.o.requires
src/CMakeFiles/lib.dir/requires: src/CMakeFiles/lib.dir/mod/f90getopt.f90.o.requires

.PHONY : src/CMakeFiles/lib.dir/requires

src/CMakeFiles/lib.dir/clean:
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src && $(CMAKE_COMMAND) -P CMakeFiles/lib.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/lib.dir/clean

src/CMakeFiles/lib.dir/depend:
	cd /home/marcus/workspace/programa_com_malhas/RMBDQD/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marcus/workspace/programa_com_malhas/RMBDQD /home/marcus/workspace/programa_com_malhas/RMBDQD/src /home/marcus/workspace/programa_com_malhas/RMBDQD/build /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src /home/marcus/workspace/programa_com_malhas/RMBDQD/build/src/CMakeFiles/lib.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/lib.dir/depend
