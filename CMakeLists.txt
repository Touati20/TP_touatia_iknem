## Indique quelle version de cmake minimale est utilisée.
cmake_minimum_required(VERSION 3.16.3)

## Specifie le nom du projet et le langage de programmation
project(TP_Particules CXX)

## Standard C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Mode Release par defaut (optimisations max)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

## Flags selon le mode :
##  - Release : ultra-rapide pour la prod
##  - Profile : avec -pg pour gprof (un peu plus lent mais profilable)
##  - Debug   : avec symboles pour gdb
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -Wall -Wextra")

## Mode Profile : on force les flags car CMake ne reconnait pas ce mode par defaut
if(CMAKE_BUILD_TYPE STREQUAL "Profile")
    add_compile_options(-O2 -pg -Wall -Wextra)
    add_link_options(-pg)
    message(STATUS "Mode PROFILE active : compilation avec -pg pour gprof")
endif()

# GoogleTest
include(FetchContent)
FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.17.0
)
FetchContent_MakeAvailable(googletest)

enable_testing()

if(APPLE)
    foreach(_gtest_target IN ITEMS gtest gtest_main gmock gmock_main)
        if(TARGET ${_gtest_target})
            set_property(
                TARGET ${_gtest_target}
                PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ""
            )
        endif()
    endforeach()
endif()

## Localisation des entetes
include_directories(${CMAKE_SOURCE_DIR}/include)

## Sous-repertoires
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(demo)
add_subdirectory(doc)