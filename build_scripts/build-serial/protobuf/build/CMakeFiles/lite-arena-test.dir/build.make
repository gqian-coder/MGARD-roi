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
CMAKE_SOURCE_DIR = /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/cmake

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build

# Include any dependencies generated for this target.
include CMakeFiles/lite-arena-test.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/lite-arena-test.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/lite-arena-test.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/lite-arena-test.dir/flags.make

/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc: protoc-3.19.4.0
/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.proto
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc"
	./protoc /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.proto --proto_path=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --cpp_out=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --experimental_allow_proto3_optional

/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc: protoc-3.19.4.0
/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.proto
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc"
	./protoc /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.proto --proto_path=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --cpp_out=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --experimental_allow_proto3_optional

/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc: protoc-3.19.4.0
/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.proto
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc"
	./protoc /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.proto --proto_path=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --cpp_out=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --experimental_allow_proto3_optional

/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc: protoc-3.19.4.0
/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.proto
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Generating /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc"
	./protoc /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.proto --proto_path=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --cpp_out=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src --experimental_allow_proto3_optional

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.s

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o: CMakeFiles/lite-arena-test.dir/flags.make
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc
CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o: CMakeFiles/lite-arena-test.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o -MF CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o.d -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o -c /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc > CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.i

CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc -o CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.s

# Object files for target lite-arena-test
lite__arena__test_OBJECTS = \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o" \
"CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o"

# External object files for target lite-arena-test
lite__arena__test_EXTERNAL_OBJECTS =

lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/lite_arena_unittest.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/arena_test_util.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_test_util.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/test_util_lite.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc.o
lite-arena-test: CMakeFiles/lite-arena-test.dir/build.make
lite-arena-test: libprotobuf-lite.so.3.19.4.0
lite-arena-test: libgmock_main.a
lite-arena-test: libgmock.a
lite-arena-test: CMakeFiles/lite-arena-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Linking CXX executable lite-arena-test"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lite-arena-test.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/lite-arena-test.dir/build: lite-arena-test
.PHONY : CMakeFiles/lite-arena-test.dir/build

CMakeFiles/lite-arena-test.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/lite-arena-test.dir/cmake_clean.cmake
.PHONY : CMakeFiles/lite-arena-test.dir/clean

CMakeFiles/lite-arena-test.dir/depend: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/map_lite_unittest.pb.cc
CMakeFiles/lite-arena-test.dir/depend: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_lite.pb.cc
CMakeFiles/lite-arena-test.dir/depend: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_import_public_lite.pb.cc
CMakeFiles/lite-arena-test.dir/depend: /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/src/google/protobuf/unittest_lite.pb.cc
	cd /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/cmake /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/src/cmake /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build /home/qian/Software/MGARD-adaptive/build_scripts/build-serial/protobuf/build/CMakeFiles/lite-arena-test.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/lite-arena-test.dir/depend

