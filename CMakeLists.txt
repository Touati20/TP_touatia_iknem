## Indique qu'elle version de cmake minimale est utilisée.
cmake_minimum_required (VERSION 3.16.3)

## Specifie le nom du projet et le langage de programmation du projet
project(Lab5 CXX)

## Force le standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Mode Release par défaut (= optimisations max si CMAKE_BUILD_TYPE non spécifié)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

## Options de compilation : -O3 (optims max) -march=native (CPU spécifique) -DNDEBUG (pas d'asserts)
set(CMAKE_CXX_FLAGS "-O3 -march=native -DNDEBUG -Wall -Wextra")

# === GoogleTest ===
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

## Définir la localisation des entêtes
include_directories(${CMAKE_SOURCE_DIR}/include)

## Sous-répertoires
add_subdirectory(src)
add_subdirectory(test)
add_subdirectory(demo)