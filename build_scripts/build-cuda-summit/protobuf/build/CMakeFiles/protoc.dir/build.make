# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/cmake

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build

# Include any dependencies generated for this target.
include CMakeFiles/protoc.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/protoc.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/protoc.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/protoc.dir/flags.make

CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o: CMakeFiles/protoc.dir/flags.make
CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o: /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc
CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o: CMakeFiles/protoc.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o -MF CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o.d -o CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o -c /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc

CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc > CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.i

CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc -o CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.s

# Object files for target protoc
protoc_OBJECTS = \
"CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o"

# External object files for target protoc
protoc_EXTERNAL_OBJECTS =

protoc-3.19.4.0: CMakeFiles/protoc.dir/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/src/google/protobuf/compiler/main.cc.o
protoc-3.19.4.0: CMakeFiles/protoc.dir/build.make
protoc-3.19.4.0: libprotoc.so.3.19.4.0
protoc-3.19.4.0: libprotobuf.so.3.19.4.0
protoc-3.19.4.0: /usr/lib/x86_64-linux-gnu/libz.so
protoc-3.19.4.0: CMakeFiles/protoc.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable protoc"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/protoc.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_executable protoc-3.19.4.0 protoc

protoc: protoc-3.19.4.0

# Rule to build all files generated by this target.
CMakeFiles/protoc.dir/build: protoc
.PHONY : CMakeFiles/protoc.dir/build

CMakeFiles/protoc.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/protoc.dir/cmake_clean.cmake
.PHONY : CMakeFiles/protoc.dir/clean

CMakeFiles/protoc.dir/depend:
	cd /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/cmake /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/src/cmake /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build /home/qian/Software/MGARD/build_scripts/build-cuda-summit/protobuf/build/CMakeFiles/protoc.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/protoc.dir/depend

