# Créer un executable à partir de fichiers sources


add_executable(run_tests 
    test_particules.cxx
    test_vecteur.cxx
    test_univers.cxx
)


# Fait le lien avec la bibliothèque créée précedemment
# Les entêtes sont automatiquement détectés.
target_link_libraries(run_tests GTest::gtest_main particules)

include(GoogleTest)
gtest_discover_tests(run_tests)


