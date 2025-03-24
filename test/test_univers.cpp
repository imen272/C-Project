#include "Univers.h"
#include <gtest/gtest.h>
using namespace std;

// Test de l'initialisation de l'univers
TEST(UniversTest, Initialisation) {
    vector<double> ld = {10.0, 10.0}; // Limites de l'univers
    Univers univers(27, 2, ld, 1.0); // 27 particules, dimension 2, rcut = 1.0
    EXPECT_EQ(univers.getNb(), 27); // Vérifie le nombre de particules
}

// Test du potentiel de Lennard-Jones
TEST(UniversTest, PotentielLennardJones) {
    vector<double> ld = {10.0, 10.0}; // Limites de l'univers
    Univers univers(27, 2, ld, 1.0); // 27 particules, dimension 2, rcut = 1.0
    double sigma = 1.0;
    double epsilon = 1.0;
    double potentiel = univers.potentielLennardJones(sigma, epsilon);
    EXPECT_NEAR(potentiel, 0.0, 1e-6); // Vérifie que le potentiel est nul à l'initialisation
}




