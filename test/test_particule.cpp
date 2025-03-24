#include <cmath>
#include "Particule.h"
#include <gtest/gtest.h>


TEST(ParticuleTest, ForceInteractionAtSigma) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(1, 0, 0); // distance de 1 entre les deux
    Particule p1(1, position1); //id 1
    Particule p2(2, position2); //id 2

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;

    // Calcul de la force
    Vecteur force = p1.forceInteraction(p2, sigma, epsilon);

    // Valeur attendue
    EXPECT_NEAR(force[0], 24 * epsilon, 1e-6); // Force en x
    EXPECT_NEAR(force[1], 0.0, 1e-6); // Force en y
    EXPECT_NEAR(force[2], 0.0, 1e-6); // Force en z
}

TEST(ParticuleTest, ForceInteractionAtMinimum) {
    Vecteur position1(0, 0, 0);
    double r_min = pow(2, 1.0 / 6.0); 
    Vecteur position2(r_min, 0, 0);   // Distance de r_min entre les deux particules
    Particule p1(1, position1);      
    Particule p2(2, position2);  

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;

    // Calcul de la force
    Vecteur force = p1.forceInteraction(p2, sigma, epsilon);

    // Valeur attendue : la force doit être nulle à r = 2^(1/6) * sigma
    EXPECT_NEAR(force[0], 0.0, 1e-6); // Force en x
    EXPECT_NEAR(force[1], 0.0, 1e-6); // Force en y
    EXPECT_NEAR(force[2], 0.0, 1e-6); // Force en z
}

TEST(ParticuleTest, ForceInteractionAt2Sigma) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(2, 0, 0); // Distance de 2 entre les deux particules
    Particule p1(1, position1); 
    Particule p2(2, position2);

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;

    // Calcul de la force
    Vecteur force = p1.forceInteraction(p2, sigma, epsilon);

    // Valeur attendue : 
    double expectedForceX = -0.3633 * epsilon;
    EXPECT_NEAR(force[0], expectedForceX, 1e-6); 
    EXPECT_NEAR(force[1], 0.0, 1e-6);           
    EXPECT_NEAR(force[2], 0.0, 1e-6);          
}

TEST(ParticuleTest, ForceInteractionGravitationelleAtDistance1) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(1, 0, 0); // Distance de 1 entre les deux particules
    Particule p1(1, 1.0, position1); // id = 1, masse = 1.0
    Particule p2(2, 1.0, position2); // id = 2, masse = 1.0

    // Calcul de la force
    Vecteur force = p1.forceInteractionGravitationelle(p2);

    // Valeur attendue
    const double G = 6.674e-11;
    EXPECT_NEAR(force[0], G, 1e-6); // Force en x
    EXPECT_NEAR(force[1], 0.0, 1e-6); // Force en y
    EXPECT_NEAR(force[2], 0.0, 1e-6); // Force en z
}

TEST(ParticuleTest, ForceInteractionGravitationelleAtDistance2) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(2, 0, 0); // Distance de 2 entre les deux particules
    Particule p1(1, 2.0, position1); // id = 1, masse = 2.0
    Particule p2(2, 3.0, position2); // id = 2, masse = 3.0

    // Calcul de la force
    Vecteur force = p1.forceInteractionGravitationelle(p2);

    // Valeur attendue
    const double G = 6.674e-11;
    double expectedForceX = 1.5 * G; // 1.5 * G
    EXPECT_NEAR(force[0], expectedForceX, 1e-6); // Force en x
    EXPECT_NEAR(force[1], 0.0, 1e-6); // Force en y
    EXPECT_NEAR(force[2], 0.0, 1e-6); // Force en z
}

TEST(ParticuleTest, ForceInteractionGravitationelleSelf) {
    Vecteur position1(0, 0, 0);
    Particule p1(1, 1.0, position1); // id = 1, masse = 1.0

    // Calcul de la force
    Vecteur force = p1.forceInteractionGravitationelle(p1);

    // Valeur attendue : la force doit être nulle
    EXPECT_NEAR(force[0], 0.0, 1e-6); // Force en x
    EXPECT_NEAR(force[1], 0.0, 1e-6); // Force en y
    EXPECT_NEAR(force[2], 0.0, 1e-6); // Force en z
}

TEST(ParticuleTest, PotentielInteractionAtSigma) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(1, 0, 0); // Distance de 1 entre les deux particules
    Particule p1(1, position1); // id = 1
    Particule p2(2, position2); // id = 2

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5; // rcut > sigma

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue : le potentiel doit être nul à r = sigma
    EXPECT_NEAR(potentiel, 0.0, 1e-6);
}

TEST(ParticuleTest, PotentielInteractionAtMinimum) {
    Vecteur position1(0, 0, 0);
    double r_min = pow(2, 1.0 / 6.0); // Distance minimale
    Vecteur position2(r_min, 0, 0);   // Distance de r_min entre les deux particules
    Particule p1(1, position1);      // id = 1
    Particule p2(2, position2);      // id = 2

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5; // rcut > r_min

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue : le potentiel doit être -epsilon à r = 2^(1/6) * sigma
    EXPECT_NEAR(potentiel, -epsilon, 1e-6);
}

TEST(ParticuleTest, PotentielInteractionAt2Sigma) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(2, 0, 0); // Distance de 2 entre les deux particules
    Particule p1(1, position1); // id = 1
    Particule p2(2, position2); // id = 2

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5; // rcut > 2 * sigma

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue : le potentiel doit être approximativement -0.0615234 * epsilon
    double expectedPotentiel = -0.0615234 * epsilon;
    EXPECT_NEAR(potentiel, expectedPotentiel, 1e-6);
}

TEST(ParticuleTest, PotentielInteractionBeyondRcut) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(3, 0, 0); // Distance de 3 entre les deux particules
    Particule p1(1, position1); // id = 1
    Particule p2(2, position2); // id = 2

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5; // rcut < 3

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue : le potentiel doit être nul à r > rcut
    EXPECT_NEAR(potentiel, 0.0, 1e-6);
}

TEST(ParticuleTest, PotentielInteractionAtZeroDistance) {
    Vecteur position1(0, 0, 0);
    Particule p1(1, position1); // id = 1
    Particule p2(2, position1); // id = 2, même position que p1

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5;

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue : le potentiel doit être nul à r = 0
    EXPECT_NEAR(potentiel, 0.0, 1e-6);
}
