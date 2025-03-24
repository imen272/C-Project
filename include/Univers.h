#ifndef UNIVERS_H
#define UNIVERS_H
#include "Particule.h"
#include <vector>
using namespace std;


class Cellule {
    public:
    vector<Particule> particules;
    vector<int> voisins;
    void ajouterParticule(const Particule& p) {
        particules.push_back(p);
    };
    int getSize() const{
        return particules.size();
    };
};


class Univers {
    private:
        vector<Cellule> grille; //1D
        int dimension,nb;
        double rcut;
        vector<double> ld;
        vector<int> ncd;

    public:
        Univers(int nb,int dimension, vector<double> ld, double rcut);
        Univers(int nb1, int nb2, int dimension, vector<double> ld, double rcut, double sigma);
        int getNb() const;
        int getIndice(const std::vector<int>& coords) const;
        void forceParticule(double sigma, double epsilon);
        double potentielLennardJones(double sigma, double epsilon);
        void etatUnivers(double dt,double sigma, double epsilon,double t_end);
        void initialiserCellulesVoisines();
        Vecteur getCentreCellule(int cellIndex);
        void mettreAJourPositionDansLimites(Particule& particule);
        void miseAJourCellules();
        void sauvegarderVTK(const string& nomFichier) const;
};

#endif
