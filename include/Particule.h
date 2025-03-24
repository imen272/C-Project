#include <string>
using namespace std;
#ifndef PARTICULE_H
#define PARTICULE_H
#include "Vecteur.h"


class Particule{
    private:
        Vecteur position;
        Vecteur vitesse;
        double masse;
        double identifiant;
        string type;
        Vecteur force;
        Vecteur oldForce;

    public:
        Particule();
        Particule(double id);
        Particule(double id,Vecteur& p);
        Particule(double id,double m,Vecteur& p);
        //Particule(double id,const Vecteur v,const Vecteur f, const Vecteur p, const string t, double m);
        Particule(double id,Vecteur p, Vecteur v,string t, double m);      
        void afficher();
        void updateVitesse(double& dt);
        void updatePosition(double dt);
        void setForce(Vecteur& f);
        void setVitesse(const Vecteur& v);
        void setPosition(Vecteur& p);
        double getIdentifiant();
        double getMasse() const;
        string getType() const;
        const Vecteur& getPosition() const;
        const Vecteur& getVitesse() const;
        Vecteur getOldForce() const;
        Vecteur getForce() const;
        Vecteur forceInteractionGravitationelle(Particule& p);
        Vecteur forceInteraction(Particule& p,double sigma, double epsilon);
        double potentielInteraction(Particule& p,double rcut,double sigma, double epsilon);
};

#endif
