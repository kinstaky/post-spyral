# Require the cmake minimum version.
cmake_minimum_required(VERSION 3.10)


# Set project properties.
project(
	a1975
	VERSION 0.1.0
	DESCRIPTION "A project to analysis data from a1975 experiment."
	LANGUAGES CXX
)


# Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if (EXISTS "${LOC_PATH}")
	message(FATAL_ERROR "You cannot build in a source directory (or any directory with"
						"CmakeLists.txt file). Please make a build subdirectory. Feel free"
						"to remove CMakeCache.txt and CMakeFiles.")
endif()

# Set the build type
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()


include(FetchContent)

# fetch catima
FetchContent_Declare(
	catima
	GIT_REPOSITORY https://github.com/kinstaky/catima.git
	GIT_TAG v1.90
	GIT_SHALLOW TRUE
	GIT_SHALLOW_DEPTH 1
)
# make available
FetchContent_MakeAvailable(catima)

# fetch toml11
FetchContent_Declare(
  toml11
  GIT_REPOSITORY https://github.com/ToruNiina/toml11.git
  GIT_TAG        v4.4.0
  GIT_SHALLOW TRUE
  GIT_SHALLOW_DEPTH 1
)
FetchContent_MakeAvailable(toml11)

# fetch cxxopts
FetchContent_Declare(
	cxxopts
	GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
	GIT_TAG v3.3.1
	GIT_SHALLOW TRUE
	GIT_SHALLOW_DEPTH 1
)
FetchContent_MakeAvailable(cxxopts)



# Only do these if this is the main project, and not if it is included
# through add_subdirectory
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

	# Let's ensure -std=c++xx instead of -std=g++xx
	set(CMAKE_CXX_EXTENSIONS OFF)

	# Set the default cxx standard
	set(CMAKE_CXX_STANDARD 17)

	# Set default compile options
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wextra")

	# Let's nicely support folders in IDE's
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)

	# Enable testing if it's in the main app.
	include(CTest)
endif()

# find root
find_package(ROOT 6.36 REQUIRED COMPONENTS Core RIO Tree Hist Graf)


option(BUILD_TESTING "Build test" OFF)
option(USE_HARDCODE_DATA "Use hardcode AMDC2020 atomic mass data" OFF)

# option(BUILD_ALL "Builds all parts" ON)
# option(BUILD_ANALYSIS "Builds analysis part" OFF)
# option(BUILD_NUMBERIC_OPTIMIZE "Builds numberic optmization part" OFF)
# option(BUILD_SIMULATION "Builds simulation part" OFF)

# if (BUILD_ALL)
# 	set(BUILD_ANALYSIS ON)
# 	set(BUILD_NUMBERIC_OPTIMIZE ON)
# 	set(BUILD_SIMULATION ON)
# endif()

function(show_config_options)
	message("")
	message("Build Properties")
	message("  c++ compiler.................: " ${CMAKE_CXX_COMPILER})
	message("  c++ stardard.................: " ${CMAKE_CXX_STANDARD})
	message("  cmake........................: " ${CMAKE_VERSION})
	message("  cmake build type.............: " ${CMAKE_BUILD_TYPE})
	message("")
	message("Build Options")
	message("  build testing................: " ${BUILD_TESTING})
	message("User-defined Flags")
	message("  USE_HARDCODE_DATA............: " ${USE_HARDCODE_DATA})
endfunction()

show_config_options()

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	add_compile_definitions(DEBUG)
endif()

if (USE_HARDCODE_DATA)
	add_compile_definitions(USE_HARDCODE_DATA)
endif()

add_subdirectory(src)
add_subdirectory(bin)

if (BUILD_TESTING)
	add_subdirectory(test)
endif()


file(
	COPY ${PROJECT_SOURCE_DIR}/assets/
	DESTINATION ${PROJECT_BINARY_DIR}/assets/
)