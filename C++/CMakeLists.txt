cmake_minimum_required(VERSION 3.10)

## set the project name (meaningless) and programming language
project(run_tumble_on_sphere LANGUAGES CXX)

## cmake commandline options (options must appear before source directory): cmake -DBUILD_DEBUG=ON ../src
option(ENABLE_OPENMP "Turn on parallelization" OFF)
## at max only one of these options should be true; TODO: throw error if that is not the case for DAUs
option(BUILD_PROFIlER_KCACHEGRIND "Build profiler (kcachegrind) version" OFF)
option(BUILD_PROFIlER_GPROF "Build profiler (gprof) version" OFF)
option(BUILD_DEBUG "Build segfault_finder version for valgrind" OFF)			# TODO: Add .dbg. in front of exe



## get name of machine
cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
message(STATUS "HOSTNAME: ${hostname}")

## specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

## show compilation commands on terminal
# set(CMAKE_VERBOSE_MAKEFILE ON)

## googletest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()


## OpenMP
if(ENABLE_OPENMP)
	find_package(OpenMP REQUIRED)
	# TODO: add output of openmp version
	#~ message(STATUS ${OpenMP_VERSION})
	link_libraries(OpenMP::OpenMP_CXX)
endif()


## Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
link_libraries(Eigen3::Eigen)
add_compile_options(-DEIGEN_NO_DEBUG -fno-math-errno)			# for speedup: no range check (-DEIGEN_NO_DEBUG) and no sqrt(), floor() checks (-fno-math-errno)

# static linking in Windows
if(WIN32)
	link_libraries(-static-libstdc++)
endif(WIN32)

## compilation options / CXX_FLAGS
if(BUILD_DEBUG)
	add_compile_options(-Wall -O0 -g)						# for valgrind when searching segfaults: valgrind --leak-check=full --track-origins=yes ../run_sim_test.exe params.task1.in.cfg 2> valgrindoutput.txt
elseif(BUILD_PROFIlER_KCACHEGRIND)
	add_compile_options(-Wall -O3 -g)						# for kcachegrind profiler
elseif(BUILD_PROFIlER_GPROF)
	add_compile_options(-Wall -O3 -g -pg)					# for gprof profiler
else()
	add_compile_options(-Wall -O3)							# default: for release/production mode
endif()



## source files
set(COMMON_SOURCES
	./src/utility_functions.cpp
	./src/run_tumble_on_sphere.cpp
	./src/Params.cpp
	./src/droplet.cpp
)



if(UNIX)
	set(TARGET1 "run_sim_unix.exe")
endif(UNIX)
if(WIN32)
	set(TARGET1 "run_sim_win")
endif(WIN32)
add_executable(${TARGET1} ${COMMON_SOURCES})



## unit tests
add_executable(
	test_droplet
	./src/droplet.cpp
	./unit_tests/test_droplet.cpp
)

target_link_libraries(
	test_droplet
	GTest::gmock_main
)

include(GoogleTest)
gtest_discover_tests(test_droplet)