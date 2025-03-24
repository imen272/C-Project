#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include "Vecteur.h"
#include "debug.h"

Vecteur::Vecteur(float u , float v , float w): x(u), y(v), z(w) {};
Vecteur::Vecteur(float u , float v): x(u), y(v), z(0.0f) {};
Vecteur::Vecteur(float u): x(u), y(0.0f), z(0.0f) {};

Vecteur::Vecteur(): x(0.0f), y(0.0f), z(0.0f){};
Vecteur Vecteur::operator+(const Vecteur& v) const {
    return Vecteur(x + v.x, y + v.y, z + v.z);
}
Vecteur Vecteur::operator-(const Vecteur& v) const {
    return Vecteur(x - v.x, y - v.y, z - v.z); 
}
Vecteur Vecteur::operator*(double scalar) const { 
    return Vecteur(x * scalar, y * scalar, z * scalar); 
}

Vecteur Vecteur::operator/(double scalar) const { 
    return Vecteur(floor(x / scalar), floor(y / scalar), floor(z / scalar)); 
}

double Vecteur::norm() const {
    return sqrt(x*x + y*y + z*z); 
}
Vecteur& Vecteur::operator+=(const Vecteur& v) {
    this->x += v.x;
    this->y += v.y;
    this->z += v.z;
    return *this;
}
/*float Vecteur:: getx() const{
    return this->x;
};
float Vecteur:: gety() const{
    return this->y;
};
float Vecteur:: getz() const{
    return this->z;
};*/

// Operator[] pour acces aux elements
float Vecteur::operator[](int index) const {
    switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: throw std::out_of_range("Index out of range for Vecteur");
    }
}

float& Vecteur::operator[](int index) {
    switch (index) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default: throw std::out_of_range("Index out of range for Vecteur");
    }
}

void Vecteur:: afficher(){
  DEBUG_COUT( "Les coordonnées du vecteur sont : " <<  this->x <<',' << this->y <<',' << this->z << endl);
}
