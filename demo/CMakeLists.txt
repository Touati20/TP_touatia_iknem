# demo/CMakeLists.txt
# Crée l'exécutable de démonstration : simulation d'impact
# d'un disque de Lennard-Jones sur une membrane.

add_executable(demo demo.cpp)

# Lien avec la bibliothèque particules (Particule + Univers + vecteur)
target_link_libraries(demo PRIVATE particules)
target_compile_features(demo PRIVATE cxx_std_17)

add_custom_command(TARGET demo POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
            ${CMAKE_BINARY_DIR}/frames
)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.0)
    target_link_libraries(demo PRIVATE stdc++fs)
endif()