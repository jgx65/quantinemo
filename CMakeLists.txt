cmake_minimum_required(VERSION 3.16)

project(quantiNemo
        VERSION 2.0.2
        DESCRIPTION "Individual-based stochastic simulation for population and quantitative genetics"
        LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Language standard
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Default to a Release build when the user did not pick one.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo)
endif()

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
option(QN_STATIC "Link the executable statically" OFF)
option(QN_THREAD "Enable the multithreaded sweep (-D_THREAD)" OFF)

# ---------------------------------------------------------------------------
# Sources — globbed so this list can never drift from the Makefile's glob.
# CONFIGURE_DEPENDS re-runs the glob when files are added/removed.
# ---------------------------------------------------------------------------
file(GLOB SRCS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB HEADERS CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

add_executable(quantinemo ${SRCS} ${HEADERS})

target_include_directories(quantinemo PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_compile_definitions(quantinemo PRIVATE _SHOW_MEMORY)
if(QN_THREAD)
    target_compile_definitions(quantinemo PRIVATE _THREAD)
    find_package(Threads REQUIRED)
    target_link_libraries(quantinemo PRIVATE Threads::Threads)
endif()

# ---------------------------------------------------------------------------
# Embed the git version (mirrors the Makefile's VERSIONGIT define).
# ---------------------------------------------------------------------------
find_package(Git QUIET)
if(GIT_FOUND)
    execute_process(
        COMMAND "${GIT_EXECUTABLE}" describe --abbrev=7 --always --tags --dirty
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        OUTPUT_VARIABLE QN_GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET)
endif()
if(QN_GIT_VERSION)
    target_compile_definitions(quantinemo PRIVATE "VERSIONGIT=\"${QN_GIT_VERSION}\"")
endif()

# ---------------------------------------------------------------------------
# Static linking (opt-in; off by default so it links on macOS).
# ---------------------------------------------------------------------------
if(QN_STATIC)
    if(APPLE)
        message(WARNING "QN_STATIC is not supported on macOS; ignoring.")
    else()
        target_link_options(quantinemo PRIVATE -static)
    endif()
endif()

# Put the binary in bin/ to match the Makefile layout.
set_target_properties(quantinemo PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
