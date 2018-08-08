# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Default target executed when no arguments are given to make.
default_target: all

.PHONY : default_target

# Allow only one "make -f Makefile2" at a time, but pass parallelism.
.NOTPARALLEL:


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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/s23410fw/Documents/software/coinfinder/coinfinder-plus/coinfinder-plus

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/s23410fw/Documents/software/coinfinder/coinfinder-plus/coinfinder-plus

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/local/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache

.PHONY : rebuild_cache/fast

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake cache editor..."
	/usr/local/bin/ccmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache

.PHONY : edit_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/s23410fw/Documents/software/coinfinder/coinfinder-plus/coinfinder-plus/CMakeFiles /Users/s23410fw/Documents/software/coinfinder/coinfinder-plus/coinfinder-plus/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /Users/s23410fw/Documents/software/coinfinder/coinfinder-plus/coinfinder-plus/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean

.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named coinfinder

# Build rule for target.
coinfinder: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 coinfinder
.PHONY : coinfinder

# fast build rule for target.
coinfinder/fast:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/build
.PHONY : coinfinder/fast

binomial_test.o: binomial_test.cpp.o

.PHONY : binomial_test.o

# target to build an object file
binomial_test.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/binomial_test.cpp.o
.PHONY : binomial_test.cpp.o

binomial_test.i: binomial_test.cpp.i

.PHONY : binomial_test.i

# target to preprocess a source file
binomial_test.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/binomial_test.cpp.i
.PHONY : binomial_test.cpp.i

binomial_test.s: binomial_test.cpp.s

.PHONY : binomial_test.s

# target to generate assembly for a file
binomial_test.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/binomial_test.cpp.s
.PHONY : binomial_test.cpp.s

coincidence.o: coincidence.cpp.o

.PHONY : coincidence.o

# target to build an object file
coincidence.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/coincidence.cpp.o
.PHONY : coincidence.cpp.o

coincidence.i: coincidence.cpp.i

.PHONY : coincidence.i

# target to preprocess a source file
coincidence.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/coincidence.cpp.i
.PHONY : coincidence.cpp.i

coincidence.s: coincidence.cpp.s

.PHONY : coincidence.s

# target to generate assembly for a file
coincidence.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/coincidence.cpp.s
.PHONY : coincidence.cpp.s

connectivity.o: connectivity.cpp.o

.PHONY : connectivity.o

# target to build an object file
connectivity.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/connectivity.cpp.o
.PHONY : connectivity.cpp.o

connectivity.i: connectivity.cpp.i

.PHONY : connectivity.i

# target to preprocess a source file
connectivity.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/connectivity.cpp.i
.PHONY : connectivity.cpp.i

connectivity.s: connectivity.cpp.s

.PHONY : connectivity.s

# target to generate assembly for a file
connectivity.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/connectivity.cpp.s
.PHONY : connectivity.cpp.s

dataset.o: dataset.cpp.o

.PHONY : dataset.o

# target to build an object file
dataset.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/dataset.cpp.o
.PHONY : dataset.cpp.o

dataset.i: dataset.cpp.i

.PHONY : dataset.i

# target to preprocess a source file
dataset.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/dataset.cpp.i
.PHONY : dataset.cpp.i

dataset.s: dataset.cpp.s

.PHONY : dataset.s

# target to generate assembly for a file
dataset.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/dataset.cpp.s
.PHONY : dataset.cpp.s

elements.o: elements.cpp.o

.PHONY : elements.o

# target to build an object file
elements.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/elements.cpp.o
.PHONY : elements.cpp.o

elements.i: elements.cpp.i

.PHONY : elements.i

# target to preprocess a source file
elements.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/elements.cpp.i
.PHONY : elements.cpp.i

elements.s: elements.cpp.s

.PHONY : elements.s

# target to generate assembly for a file
elements.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/elements.cpp.s
.PHONY : elements.cpp.s

id_lookup.o: id_lookup.cpp.o

.PHONY : id_lookup.o

# target to build an object file
id_lookup.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/id_lookup.cpp.o
.PHONY : id_lookup.cpp.o

id_lookup.i: id_lookup.cpp.i

.PHONY : id_lookup.i

# target to preprocess a source file
id_lookup.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/id_lookup.cpp.i
.PHONY : id_lookup.cpp.i

id_lookup.s: id_lookup.cpp.s

.PHONY : id_lookup.s

# target to generate assembly for a file
id_lookup.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/id_lookup.cpp.s
.PHONY : id_lookup.cpp.s

main.o: main.cpp.o

.PHONY : main.o

# target to build an object file
main.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/main.cpp.o
.PHONY : main.cpp.o

main.i: main.cpp.i

.PHONY : main.i

# target to preprocess a source file
main.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/main.cpp.i
.PHONY : main.cpp.i

main.s: main.cpp.s

.PHONY : main.s

# target to generate assembly for a file
main.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/main.cpp.s
.PHONY : main.cpp.s

parameters.o: parameters.cpp.o

.PHONY : parameters.o

# target to build an object file
parameters.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/parameters.cpp.o
.PHONY : parameters.cpp.o

parameters.i: parameters.cpp.i

.PHONY : parameters.i

# target to preprocess a source file
parameters.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/parameters.cpp.i
.PHONY : parameters.cpp.i

parameters.s: parameters.cpp.s

.PHONY : parameters.s

# target to generate assembly for a file
parameters.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/parameters.cpp.s
.PHONY : parameters.cpp.s

significance.o: significance.cpp.o

.PHONY : significance.o

# target to build an object file
significance.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/significance.cpp.o
.PHONY : significance.cpp.o

significance.i: significance.cpp.i

.PHONY : significance.i

# target to preprocess a source file
significance.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/significance.cpp.i
.PHONY : significance.cpp.i

significance.s: significance.cpp.s

.PHONY : significance.s

# target to generate assembly for a file
significance.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/significance.cpp.s
.PHONY : significance.cpp.s

test_cases.o: test_cases.cpp.o

.PHONY : test_cases.o

# target to build an object file
test_cases.cpp.o:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/test_cases.cpp.o
.PHONY : test_cases.cpp.o

test_cases.i: test_cases.cpp.i

.PHONY : test_cases.i

# target to preprocess a source file
test_cases.cpp.i:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/test_cases.cpp.i
.PHONY : test_cases.cpp.i

test_cases.s: test_cases.cpp.s

.PHONY : test_cases.s

# target to generate assembly for a file
test_cases.cpp.s:
	$(MAKE) -f CMakeFiles/coinfinder.dir/build.make CMakeFiles/coinfinder.dir/test_cases.cpp.s
.PHONY : test_cases.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... rebuild_cache"
	@echo "... edit_cache"
	@echo "... coinfinder"
	@echo "... binomial_test.o"
	@echo "... binomial_test.i"
	@echo "... binomial_test.s"
	@echo "... coincidence.o"
	@echo "... coincidence.i"
	@echo "... coincidence.s"
	@echo "... connectivity.o"
	@echo "... connectivity.i"
	@echo "... connectivity.s"
	@echo "... dataset.o"
	@echo "... dataset.i"
	@echo "... dataset.s"
	@echo "... elements.o"
	@echo "... elements.i"
	@echo "... elements.s"
	@echo "... id_lookup.o"
	@echo "... id_lookup.i"
	@echo "... id_lookup.s"
	@echo "... main.o"
	@echo "... main.i"
	@echo "... main.s"
	@echo "... parameters.o"
	@echo "... parameters.i"
	@echo "... parameters.s"
	@echo "... significance.o"
	@echo "... significance.i"
	@echo "... significance.s"
	@echo "... test_cases.o"
	@echo "... test_cases.i"
	@echo "... test_cases.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

