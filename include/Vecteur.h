#ifndef VECTEUR_H
#define VECTEUR_H
#include <vector>
using namespace std;

class Vecteur {
    public:
        float x;
        float y;
        float z;
        Vecteur(float u, float v , float w);
        Vecteur(float u , float v);
        Vecteur(float u);
        Vecteur();
        Vecteur operator+(const Vecteur& v) const;
        Vecteur operator-(const Vecteur& v) const;
        Vecteur operator*(double scalar) const;
        Vecteur operator/(double scalar) const;
        Vecteur& operator+=(const Vecteur& v);
        bool operator==(const Vecteur& other) const;
        double norm() const ;
        float operator[](int index) const;
        float& operator[](int index) ;
        void afficher();
        //TODO: privacy pour Vecteur
};

#endif
