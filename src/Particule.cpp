#include "Particule.h"
#include "Vecteur.h"
#include <string>
#include <random>
#include <iostream>
#include <vector>
#include <set>
#include <array>
#include <list>
#include <deque>
#include <chrono>
#include "debug.h"

Particule::Particule() : identifiant(0), masse(1), type("default"),position(Vecteur()), vitesse(Vecteur()), force(Vecteur()),oldForce(Vecteur()) {};

Particule::Particule(double id,Vecteur& p):position(p), vitesse(Vecteur(1,0,0)), force(Vecteur()),oldForce({Vecteur()}){
    this->identifiant=id;
    masse=1;
    type="default";
}; //masse=1 par defaut!!!!

Particule::Particule(double id,double m,Vecteur& p):position(p), masse(m), vitesse(Vecteur(1,0,0)), force(Vecteur()),oldForce({Vecteur()}){
    this->identifiant=id;
    type="default";
};

Particule::Particule(double id):position(Vecteur()), vitesse(Vecteur()), force(Vecteur()),oldForce({Vecteur()}){
    this->identifiant=id;
    masse=1;
    type="default";
};

Particule::Particule(double id,Vecteur p, Vecteur v,string t, double m)
: identifiant(id),type(t), masse(m), position(p), vitesse(v), force(Vecteur()), oldForce(Vecteur()) {};

/*Particule::Particule(double id,const Vecteur v,const Vecteur f, const Vecteur p, const string t, double m) 
: identifiant(id), masse(m), type(t),position(p), force(f), vitesse(v), oldForce(Vecteur()) {}; //??????
*/
void Particule::updatePosition(double dt){
    /*position.x=position.x+ dt*(vitesse.x+ 0.5/masse *force.x*dt);
    position.y=position.y+ dt*(vitesse.y+ 0.5/masse *force.y*dt);
    position.z=position.z+ dt*(vitesse.z+ 0.5/masse *force.z*dt);*/
    position=position+vitesse*dt+force*(0.5*dt*dt/masse);
       
};

void Particule::updateVitesse(double& dt){
    /*vitesse.x=vitesse.x+ dt * 0.5/masse * (force.x+oldForce.x);
    vitesse.y=vitesse.y+ dt * 0.5/masse * (force.y+oldForce.y);
    vitesse.z=vitesse.z+ dt * 0.5/masse * (force.z+oldForce.z);*/
    vitesse=vitesse+(force+oldForce)*(0.5*dt/masse);
};

void Particule::setPosition(Vecteur& p){
    position=p;
};

void Particule::setVitesse(const Vecteur& v){
    vitesse=v;

};

void Particule::setForce(Vecteur& f){
    oldForce=force;
    force=f;
    //const to not allow modification? optional
};
double Particule::getIdentifiant(){
    return this->identifiant;
};

double Particule::getMasse() const{
    return this->masse;
};

string Particule::getType() const{
    return this->type;
};

const Vecteur& Particule::getPosition() const{
    return position;
};

const Vecteur& Particule::getVitesse() const{
    return vitesse;
};

Vecteur Particule::getOldForce() const{
    return this->oldForce;
};

Vecteur Particule::getForce() const{
    return this->force;
};

void Particule::afficher() {
  DEBUG_COUT( "Particule " << this->identifiant
        << " positionée au [" << this->position.x << "," << this->position.y << "," << this->position.z << "]"
	       << " de vitesse [" << this->vitesse.x << "," << this->vitesse.y << "," << this->vitesse.z << "]" <<endl);
};

Vecteur Particule::forceInteractionGravitationelle(Particule& p){
    const double G = 6.674e-11; // Constante gravitationnelle
    Vecteur r=this->position-p.getPosition();
    double d=r.norm();
    
    if (d == 0) {
        return Vecteur();
    }

    double module= G * this->masse * p.getMasse() / (d*d*d);
    return r*module;
};


Vecteur Particule::forceInteraction(Particule& p,double sigma, double epsilon){
    Vecteur r=this->position-p.getPosition();
    double d=r.norm();
    if (d==0){
        return Vecteur();
    }
    double u = pow(sigma/d,6);
    double module= 24* epsilon /(d*d) * u* (1-2*u);
    return r*module;


};//en remplacement de la force gravitationnelle ou en addition


double Particule::potentielInteraction(Particule& p,double rcut,double sigma, double epsilon){
    Vecteur r=this->position-p.getPosition();
    double d=r.norm();
    if (d<=rcut && d!= 0){
        double u = pow(sigma/d,6);
        return 4* epsilon* u* (u-1);
    } 
    return 0.0;
    
};


/*void forceParticule(vector<Particule>& listParticule){
    //first method: non optimisée
    for (auto& p:listParticule){
        Vecteur force;
        for (auto& p2:listParticule){
            if (&p != &p2) {
                force+=p.forceInteractionGravitationelle(p2);
                }
        } 
        p.setForce(force);
    }

   //New Method: Division du temps de calcul par 2
    for (auto& p : listParticule) {
        Vecteur reset(0, 0, 0);
        p.setForce(reset);
    }

    for (size_t i = 0; i < listParticule.size(); ++i) {
        for (size_t j = i + 1; j < listParticule.size(); ++j) {
            Vecteur f = listParticule[i].forceInteractionGravitationelle(listParticule[j]);
            
            Vecteur newForce1 = listParticule[i].getForce() + f;
            listParticule[i].setForce(newForce1);
            
            Vecteur newForce2 = listParticule[j].getForce() - f;
            listParticule[j].setForce(newForce2);
        }
    }
}
    

int main() {
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    
    vector<Particule> particles {
        Particule(dist(mt), Vecteur(0.0, 0.0), Vecteur(0.0, 0.0),"Soleil", 1.0),
        Particule(dist(mt), Vecteur(0.0, 1.0), Vecteur(-1.0, 0.0),"Terre", 3.0e-6),
        Particule(dist(mt), Vecteur(0.0, 5.36), Vecteur(-0.425, 0.0),"Jupiter", 9.55e-4),
        Particule(dist(mt), Vecteur(34.75, 0.0), Vecteur(0.0, 0.0296),"Halley", 1.0e-14)
    };

    // Störmer-Verlet for particles
    double t = 0.0;
    double t_end = 468.5;
    double dt= 0.015;
    forceParticule(particles);

    for (auto &p : particles) {
        p.afficher();
    };

    while (t < t_end) {
        t+= dt;
        for (auto &p : particles) {
            p.updatePosition(dt);
        };
        forceParticule(particles);
        for (auto &p : particles) {
            p.updateVitesse(dt); 
        };
        
    };

    for (auto &p : particles) {
        p.afficher();
    };

    return 0;

}*/










