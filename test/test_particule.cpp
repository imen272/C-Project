#include <cmath>
#include "Particule.h"
#include <gtest/gtest.h>
#include <chrono>
using namespace std;

TEST(ParticuleTest, ConstructeurParDefaut) {
    Particule p;
    EXPECT_EQ(p.getIdentifiant(), 0);
    EXPECT_EQ(p.getMasse(), 1);
    EXPECT_EQ(p.getType(), "default");
    EXPECT_EQ(p.getPosition(), Vecteur());
    EXPECT_EQ(p.getVitesse(), Vecteur());
}

TEST(ParticuleTest, ConstructeurAvecPosition) {
    Vecteur pos(1, 2, 3);
    Particule p(1, pos);
    EXPECT_EQ(p.getIdentifiant(), 1);
    EXPECT_EQ(p.getPosition(), pos);
    EXPECT_EQ(p.getVitesse(), Vecteur(1, 0, 0)); 
}

TEST(ParticuleTest, SettersAndGetters) {
    Particule p;
    Vecteur pos(1, 2, 3);
    Vecteur vit(0.5, 0, 0);
    Vecteur force(0, 0, 1);

    p.setPosition(pos);
    p.setVitesse(vit);
    p.setForce(force);

    EXPECT_EQ(p.getPosition(), pos);
    EXPECT_EQ(p.getVitesse(), vit);
    EXPECT_EQ(p.getForce(), force);
    EXPECT_EQ(p.getOldForce(), Vecteur()); 
}


TEST(ParticuleTest, ForceInteractionAtSigma) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(1, 0, 0); // distance de 1 entre les deux
    Particule p1(1, position1); 
    Particule p2(2, position2); 

    // Paramètres
    double sigma = 1.0;
    double epsilon = 1.0;

    // Calcul de la force
    Vecteur force = p1.forceInteraction(p2, sigma, epsilon);

    // Valeur attendue
    EXPECT_NEAR(force[0], -24 , 1e-6); 
    EXPECT_NEAR(force[1], 0.0, 1e-6); 
    EXPECT_NEAR(force[2], 0.0, 1e-6);
}

TEST(ParticuleTest, ForceInteractionAtMinimum) {
    Vecteur position1(0, 0, 0);
    double r_min = pow(2, 1.0 / 6.0); 
    Vecteur position2(r_min, 0, 0);   // Distance de r_min 
    Particule p1(1, position1);      
    Particule p2(2, position2);  

    // Paramètres de Lennard-Jones
    double sigma = 1.0;
    double epsilon = 1.0;

    // Calcul de la force
    Vecteur force = p1.forceInteraction(p2, sigma, epsilon);

    // Valeur attendue : la force doit être nulle à r = 2^(1/6) * sigma
    EXPECT_NEAR(force[0], 0.0, 1e-6); 
    EXPECT_NEAR(force[1], 0.0, 1e-6); 
    EXPECT_NEAR(force[2], 0.0, 1e-6); 
}

TEST(ParticuleTest, ForceInteractionAt2Sigma) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(2, 0, 0); // distance = 2*sigma
    Particule p1(1, position1);
    Particule p2(2, position2);

    double sigma = 1.0;
    double epsilon = 1.0;

    Vecteur force = p1.forceInteraction(p2, sigma, epsilon);
    
    // Calcul manuel:
    // u = (1/2)^6 = 0.015625
    // module = 24*1/(4)*0.015625*(1-2*0.015625) = 6*0.015625*0.96875 ≈ 0.0908203
    // force ≈ (0.0908203*2, 0, 0)
    EXPECT_NEAR(force[0],0.181641 , 1e-6);
    EXPECT_NEAR(force[1], 0.0, 1e-6);
    EXPECT_NEAR(force[2], 0.0, 1e-6);
}

TEST(ParticuleTest, ForceInteractionGravitationelleAtDistance1) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(1, 0, 0); // Distance de 1 
    Particule p1(1, 1.0, position1); //masse = 1.0
    Particule p2(2, 1.0, position2); // masse = 1.0

    // Calcul de la force
    Vecteur force = p1.forceInteractionGravitationelle(p2);

    // Valeur attendue
    const double G = 6.674e-11;
    EXPECT_NEAR(force[0], G, 1e-6); 
    EXPECT_NEAR(force[1], 0.0, 1e-6); 
    EXPECT_NEAR(force[2], 0.0, 1e-6); 
}

TEST(ParticuleTest, ForceInteractionGravitationelleAtDistance2) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(2, 0, 0); // Distance de 2 entre les deux particules
    Particule p1(1, 2.0, position1); 
    Particule p2(2, 3.0, position2);

    // Calcul de la force
    Vecteur force = p1.forceInteractionGravitationelle(p2);

    // Valeur attendue
    const double G = 6.674e-11;
    double expectedForceX = 1.5 * G; 
    EXPECT_NEAR(force[0], expectedForceX, 1e-6); 
    EXPECT_NEAR(force[1], 0.0, 1e-6); 
    EXPECT_NEAR(force[2], 0.0, 1e-6); 
}

TEST(ParticuleTest, ForceGravitationnelleNulleSiMemePosition) {
    Particule p1(1), p2(2);
    Vecteur force = p1.forceInteractionGravitationelle(p2);
    EXPECT_EQ(force, Vecteur());
}

TEST(ParticuleTest, PotentielLennardJonesCutoff) {
    Vecteur pos1(0, 0, 0);
    Vecteur pos2(10, 0, 0); // Distance > rcut
    Particule p1(1, pos1);
    Particule p2(2, pos2);
    
    double potentiel = p1.potentielInteraction(p2, 5.0, 1.0, 1.0);
    EXPECT_EQ(potentiel, 0.0);
}

TEST(ParticuleTest, PotentielLennardJonesAtSigma) {
    Vecteur pos1(0, 0, 0);
    Vecteur pos2(1, 0, 0); // distance = sigma
    Particule p1(1, pos1);
    Particule p2(2, pos2);
    
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 1.1; // > sigma

    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);
    
    // Calcul manuel:
    // u = (1/1)^6 = 1
    // potentiel = 4*1*1*(1-1) = 0
    EXPECT_NEAR(potentiel, 0.0, 1e-6);
}

TEST(ParticuleTest, PotentielInteractionBeyondRcut) {
    Vecteur position1(0, 0, 0);
    Vecteur position2(3, 0, 0); // Distance de 3 
    Particule p1(1, position1); 
    Particule p2(2, position2); 

    // Paramètres
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5; // rcut < 3

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue 
    EXPECT_NEAR(potentiel, 0.0, 1e-6);
}

TEST(ParticuleTest, PotentielInteractionAtZeroDistance) {
    Vecteur position1(0, 0, 0);
    Particule p1(1, position1); 
    Particule p2(2, position1); // même position que p1

    // Paramètres 
    double sigma = 1.0;
    double epsilon = 1.0;
    double rcut = 2.5;

    // Calcul du potentiel
    double potentiel = p1.potentielInteraction(p2, rcut, sigma, epsilon);

    // Valeur attendue 
    EXPECT_NEAR(potentiel, 0.0, 1e-6);
}


TEST(ParticuleTest, UpdateForcePreserveOldForce) {
    Particule p;
    Vecteur f1(1, 0, 0);
    Vecteur f2(0, 1, 0);
    
    p.setForce(f1);
    p.setForce(f2);
    
    EXPECT_EQ(p.getOldForce(), f1);
    EXPECT_EQ(p.getForce(), f2);
}

TEST(ParticuleTest, MasseNulle) {
    Particule p(123); 
    Vecteur v; 
    p.setForce(v);
    
    //devrait gérer la division par zéro
    p.updatePosition(1.0);
    p.updateVitesse(1.0);
    
    SUCCEED();  //Marqueur de réussite
}

TEST(ParticuleTest, PerformanceForceCalculation) {
    const int N = 1000;
    vector<Particule> particules;
    for (int i = 0; i < N; ++i) {
        particules.emplace_back(i, Vecteur(i, 0, 0));
    }
    
    auto start = chrono::high_resolution_clock::now();
    for (auto& p1 : particules) {
        for (auto& p2 : particules) {
            p1.forceInteractionGravitationelle(p2);
        }
    }
    auto end = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    EXPECT_LT(duration.count(), 100); 
}

