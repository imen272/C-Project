#include "Vecteur.h"
#include <gtest/gtest.h>

TEST(VecteurTest, Norme) {
    Vecteur v(3, 4, 0);
    EXPECT_EQ(v.norm(), 5);
}

TEST(VecteurTest, Soustraction) {
    Vecteur v1(5, 5, 5);
    Vecteur v2(2, 3, 4);
    Vecteur result = v1 - v2;
    EXPECT_EQ(result[0], 3);
    EXPECT_EQ(result[1], 2);
    EXPECT_EQ(result[2], 1);
}

TEST(VecteurTest, Addition) {
    Vecteur v1(1, 2, 3);
    Vecteur v2(4, 5, 6);
    Vecteur result = v1 + v2;
    EXPECT_EQ(result[0], 5); // Composante x
    EXPECT_EQ(result[1], 7); // Composante y
    EXPECT_EQ(result[2], 9); // Composante z
}
