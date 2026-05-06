# Generation de la documentation Doxygen
# Usage : depuis build/, lancer "make doc"

find_package(Doxygen)

if(DOXYGEN_FOUND)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/doxyfilerc
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        COMMENT "Generation de la documentation Doxygen..."
        VERBATIM
    )
    message(STATUS "Doxygen trouve : cible 'make doc' disponible")
else()
    message(WARNING "Doxygen non trouve : la cible 'doc' ne sera pas creee")
endif()
