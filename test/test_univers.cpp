#include "Univers.h"
#include <gtest/gtest.h>
using namespace std;
#include <chrono>

TEST(UniversTest, ConstructionValide) {
    vector<double> dimensions = {10.0, 10.0};
    Univers univers(100, 2, dimensions, 2.0);
    
    // Vérification via l'interface publique
    EXPECT_EQ(univers.getNb(), 100);
}

TEST(UniversTest, ConstructeurCollision) {
    vector<double> dimensions = {50.0, 50.0};
    Univers univers(1600, 6400, 2, dimensions, 1.0, 0.5);
    
    // Vérification indirecte via le potentiel
    EXPECT_NO_THROW(univers.potentielLennardJones(1.0, 1.0));
}


TEST(UniversTest, StabiliteLongTerme) {
    vector<double> dimensions = {30.0, 30.0};
    Univers univers(50, 2, dimensions, 3.0);
    
    EXPECT_NO_THROW(univers.etatUnivers(0.01, 1.0, 1.0, 1.0)); // 100 pas de temps
}

TEST(UniversTest, PerformanceCalculForces) {
    vector<double> dimensions = {20.0, 20.0};
    Univers univers(100, 2, dimensions, 2.5);
    
    auto start = chrono::high_resolution_clock::now();
    univers.forceParticule(1.0, 1.0);
    auto end = chrono::high_resolution_clock::now();
    
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    EXPECT_LT(duration.count(), 100); // Doit être rapide
}

