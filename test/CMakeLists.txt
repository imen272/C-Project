# Créer un exécutable pour chaque test
add_executable(test_vecteur test_vecteur.cpp)
add_executable(test_particule test_particule.cpp)
add_executable(test_univers test_univers.cpp)

# Lier les exécutables à la bibliothèque principale et à Google Test
target_link_libraries(test_vecteur gtest_main univers)
target_link_libraries(test_particule gtest_main univers)
target_link_libraries(test_univers gtest_main univers)

# Activer la découverte des tests
include(GoogleTest)
gtest_discover_tests(test_vecteur)
gtest_discover_tests(test_particule)
gtest_discover_tests(test_univers)
