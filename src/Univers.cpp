#include <fstream>
#include <iomanip>
#include <iostream>
#include "Particule.h"
#include "Univers.h"
#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <chrono>
#include <fstream>
#include "debug.h"

Univers::Univers(int nb,int dimension, vector<double> ld, double rcut): nb(nb),dimension(dimension), ld(ld), rcut(rcut) {
    ncd.resize(dimension);
    int totalCells = 1;

    // Calcul des ncd et Taille totale de la grille
    for (int i = 0; i < dimension; i++) {
        ncd[i] = static_cast<int>(ld[i] / rcut);
        totalCells *= ncd[i];

    }

    grille.resize(totalCells);

    // Génération aléatoire des particules
    std::random_device rd;
    std::mt19937 mt(rd());
    //std::uniform_real_distribution<double> id(0.0,10000);
    for (int i = 0; i < nb; i++) {
        std::vector<double> pos(dimension);

        for (int d = 0; d < dimension; d++) {
            uniform_real_distribution<double> dist(0.0, ld[d]); // Distribution uniforme entre 0 et ld[i]{
            pos[d] = dist(mt); // Position aléatoire dans [0, ld[d]]
        }

        Vecteur position(
            pos[0], 
            (dimension > 1) ? pos[1] : 0.0, 
            (dimension > 2) ? pos[2] : 0.0
        );

        Particule particule(i,position);
        particule.afficher();
        // Déterminer la cellule de la grille où placer la particule
        std::vector<int> cellCoords(dimension);
        for (int d = 0; d < dimension; d++) {
            cellCoords[d] = static_cast<int>(pos[d] / rcut);
        }
        
        int cellIndex = getIndice(cellCoords);
        grille[cellIndex].particules.push_back(particule);
    }
    /*for (size_t i = 0; i < grille.size(); i++) {
        DEBUG_COUT( "Cellule " << i << " contient " << grille[i].particules.size() << " particules." << endl);
    }*/
    
    initialiserCellulesVoisines();
};


//Nouveau constructeur pour la collision des deux objets
Univers::Univers(int nb1, int nb2, int dimension, vector<double> ld, double rcut, double sigma)
    : dimension(dimension), ld(ld), rcut(rcut) {
    ncd.resize(dimension);
    int totalCells = 1;

    //ncd et total number of cells
    for (int i = 0; i < dimension; i++) {
        ncd[i] = static_cast<int>(ld[i] / rcut);
        totalCells *= ncd[i];
    }

    grille.resize(totalCells);
    DEBUG_COUT( "Total cells in grille: " << totalCells << std::endl);
    // Initialiser les particules pour le carré
    double s = pow(2, 1.0 / 6.0) * sigma; // Distance entre les 2 objets!
    double distance= 0.4984;; //calculée en se basant sur les autres données du probleme
    int particlesPerSide1 = 40; // 40 particles par coté du carré
    for (int i = 0; i < particlesPerSide1; ++i) {
        for (int j = 0; j < particlesPerSide1; ++j) {
            Vecteur position( 30 + i * distance, j * distance, 0.0);
            
            //Debbug..
            /*if (position[0] < 0 || position[0] >= ld[0] || position[1] < 0 || position[1] >= ld[1]) {
                std::cerr << "Error: Particle position out of bounds: (" << position[0] << ", " << position[1] << ")" << endl;
                continue; 
            }*/

            Particule particule(i * particlesPerSide1 + j, position);
            particule.setVitesse(Vecteur(0.0, 10.0, 0.0)); // Vitesse initiale du carré
            particule.afficher();

            // Determine the cell index and add the particle to the grid
            vector<int> cellCoords(dimension);
            for (int d = 0; d < dimension; d++) {
                cellCoords[d] = static_cast<int>(position[d] / rcut);
            }
            int cellIndex = getIndice(cellCoords);
            grille[cellIndex].particules.push_back(particule);
        }
    }

    // Initialiser les particules pour le rectangle
    int particlesPerSide2X = 160; // 160 particles suivant x
    int particlesPerSide2Y = 40;  // 40 particles suivant y
    for (int i = 0; i < particlesPerSide2X; ++i) {
        for (int j = 0; j < particlesPerSide2Y; ++j) {
            Vecteur position( i * distance, 20.56 + j * distance, 0.0); //20.56 calculated too!!
            Particule particule(nb1 + i * particlesPerSide2Y + j, position);
            particule.setVitesse(Vecteur(0.0, 0.0, 0.0)); // Vitesse initiale du rectangle
            particule.afficher();

            vector<int> cellCoords(dimension);
            for (int d = 0; d < dimension; d++) {
                cellCoords[d] = static_cast<int>(position[d] / rcut);
            }
            int cellIndex = getIndice(cellCoords);
            grille[cellIndex].particules.push_back(particule);
        }
    }

    initialiserCellulesVoisines();
} //Not so generic but works



int Univers:: getNb() const{
    return this->nb;
};



int Univers::getIndice(const std::vector<int>& coords) const {
    // Convertir les coordonnées en un index 1D
    int index = 0;
    int facteur = 1;
    for (int i = 0; i < dimension; i++) {
        index += coords[i] * facteur;
        facteur *= ncd[i];
    }
    return index;
};



void Univers::initialiserCellulesVoisines() {
    // Initialise les voisins en fct de dimension
    for (int i = 0; i < grille.size(); i++) {
        vector<int> coords(dimension);
        int temp = i;

        for (int d = 0; d < dimension; d++) {
            coords[d] = temp % ncd[d];
            temp /= ncd[d];
        }
        // Affichage des coordonnées de la cellule courante
        /*DEBUG_COUT( "Cellule " << i << " a les coordonnées : ");
        for (int d = 0; d < dimension; d++) {
            DEBUG_COUT(coords[d] << ",";
        }
        DEBUG_COUT(")" << endl);*/

        // Liste des offsets possibles pour chaque dimension
        vector<int> offsets = {-1, 0, 1};

        // Créer un vecteur d'offsets pour chaque dimension
        vector<vector<int>> offsetCombinations;
        
        // Génération des combinaisons d'offsets
        // Nombre total de combinaisons
        int totalCombinations = 1;
        for (int d = 0; d < dimension; d++) {
            totalCombinations *= offsets.size();
        }

        // Génération des combinaisons d'offsets
        for (int combIdx = 0; combIdx < totalCombinations; combIdx++) {
            vector<int> currentOffsets(dimension);

            int tempComb = combIdx;
            for (int d = 0; d < dimension; d++) {
                currentOffsets[d] = offsets[tempComb % offsets.size()];
                tempComb /= offsets.size();
            }

            offsetCombinations.push_back(currentOffsets);
        }


        /*DEBUG_COUT( "Combinaisons d'offsets générées pour la cellule " << i << " : " << std::endl);
        for (const auto& combinaison : offsetCombinations) {
            for (int val : combinaison) {
                DEBUG_COUT( val << ",");
            }
            DEBUG_COUT(")" << endl);
        }*/

        // Ajoute les voisins en fonction des combinaisons d'offsets
        for (const auto& combinaison : offsetCombinations) {
            vector<int> voisinCoords = coords;

            // Applique chaque offset à la coordonnée correspondante
            for (int d = 0; d < dimension; d++) {
                voisinCoords[d] += combinaison[d];
            }

            // Vérification si les coordonnées sont valides
            bool valid = true;
            for (int j = 0; j < dimension; j++) {
                if (voisinCoords[j] < 0 || voisinCoords[j] >= ncd[j]) {
                    valid = false;
                    break;
                }
            }

            // ajoute l'indice du voisin
            if (valid) {
                int voisinIndex = getIndice(voisinCoords);
                if (voisinIndex >= 0 && voisinIndex < grille.size()) {
                    grille[i].voisins.push_back(voisinIndex);
                }
            }
        }

        //DEBUG_COUT( "Cellule " << i << " A " << grille[i].voisins.size() << " voisins." << endl);  
    }
};


Vecteur Univers::getCentreCellule(int cellIndex) {
    vector<int> coords(dimension);
    int temp = cellIndex;
    //Changement de dimensions
    for (int d = 0; d < dimension; d++) {
        coords[d] = temp % ncd[d];
        temp /= ncd[d];
    }
    //DEBUG_COUT( "Cellule " << cellIndex << " a les coordonnées (" << coords[0] << ", " << coords[1] << ")" << endl);
    //DEBUG_COUT( "ncd[0] = " << ncd[0] << ", ncd[1] = " << ncd[1] << endl);

    // Calcul des coordonnées du centre de la cellule
    Vecteur centre;
    for (int d = 0; d < dimension; d++) {
        centre[d] = (coords[d] + 0.5) * rcut; 
    }
    //DEBUG_COUT( "Centre de la cellule " << cellIndex << " : (" << centre[0] << ", " << centre[1] << ")" << endl);

    return centre;
}



double Univers::potentielLennardJones(double sigma, double epsilon){
    double potentiel=0;
    for (size_t i = 0; i < grille.size(); ++i) {
        for (auto& p:grille[i].particules){
            for (auto& indexVoisin:grille[i].voisins) {
                Vecteur v=getCentreCellule(indexVoisin)-p.getPosition();
                if (v.norm()<rcut){
                    for (auto& p2:grille[indexVoisin].particules) {
                        potentiel += p.potentielInteraction(p2,rcut,sigma,epsilon);
                        
                    }
                }
            }
        }
    }
    return potentiel;
};



void Univers::forceParticule(double sigma, double epsilon){
    for (size_t i = 0; i < grille.size(); ++i) {
        for (auto& p:grille[i].particules){
            Vecteur force;
            for (auto& indexVoisin:grille[i].voisins) {
                //DEBUG_COUT( "voisin est " << indexVoisin << endl);
                Vecteur v=getCentreCellule(indexVoisin)-p.getPosition();
                //DEBUG_COUT( "Position de la particule : (" << p.getPosition()[0] << ", " << p.getPosition()[1] << ")" << endl);
                //DEBUG_COUT( "centre du voisin est (" << v[0] << "," << v[1] << "," << v[2] << endl);
                if (v.norm()<rcut){
                    for (auto& p2:grille[indexVoisin].particules) { 
                        if (&p != &p2) {
                            force+=p.forceInteraction(p2,sigma,epsilon);
                            force+=p.forceInteractionGravitationelle(p2);
                        }
                    }
                }
            }
            p.setForce(force);
            //DEBUG_COUT( "force de la particule est (" << force[0] << "," << force[1] << "," << force[2] <<")"<< endl);
        }
    }
}; //...


void Univers::mettreAJourPositionDansLimites(Particule& particule) {
    Vecteur pos = particule.getPosition();

    for (size_t d = 0; d < dimension; d++) {  
        pos[d] = fmod(pos[d], ld[d]); 
        if (pos[d] < 0) pos[d] += ld[d]; 
    }

    particule.setPosition(pos);  
};



void Univers::miseAJourCellules() {

    vector<Particule> allParticles;
    /*for (auto& cellule : grille) {
        DEBUG_COUT( "cell has " << cellule.getSize() << endl);
   }*/

    for (auto& cellule : grille) {
        for (auto& p : cellule.particules) {
            allParticles.push_back(p); 
        }
        cellule.particules.clear();
    }
    //Reassign les particules à leurs nouvelles cellules
    for (auto& p : allParticles) {
        vector<int> newCellCoords(dimension);
        
        // Calculer les indices de la cellule pour chaque dimension
        for (int d = 0; d < dimension; d++) {
            newCellCoords[d] = static_cast<int>(floor(p.getPosition()[d] / rcut));
        }

        int cellIndex = getIndice(newCellCoords);
        grille[cellIndex].ajouterParticule(p);
    
    }
    /*for (auto& cellule : grille) {
         DEBUG_COUT( "cell has " << cellule.getSize() << endl);
    }*/
    initialiserCellulesVoisines();
};



void Univers::etatUnivers(double dt,double sigma, double epsilon,double t_end){
    double t = 0.0;
    int pas=0;
    forceParticule(sigma,epsilon);
    /*for (auto &p : listParticule) {
        p.afficher();
    };*/

    //DEBUG_COUT( "Nouvelles position et vitesse des particules :" << endl);
    while (t < t_end) {
        t+= dt;
        for (size_t i = 0; i < grille.size(); ++i) {
            for (auto &p : grille[i].particules) {
                p.updatePosition(dt);
                mettreAJourPositionDansLimites(p);
                //$$
                //p.afficher();

            }

            forceParticule(sigma,epsilon);

            for (auto &p : grille[i].particules) {
                p.updateVitesse(dt); 
                p.afficher();
            }
        }
        miseAJourCellules();

	// Sauvegarder l'état actuel dans un fichier VTK
        std::string nomFichier = "output/particules_t_" + std::to_string(pas) + ".vtu";
        sauvegarderVTK(nomFichier);
	pas+=1;
    }
};


void Univers::sauvegarderVTK(const std::string& nomFichier) const {
  string cheminComplet = string(getenv("HOME")) + "/TP_errabhii_madakcha/" +  nomFichier;
    cout << "Chemin du fichier : " << cheminComplet << endl;
    ofstream fichier(cheminComplet);
    if (!fichier.is_open()) {
        std::cerr << "Erreur : Impossible d'ouvrir le fichier " << nomFichier << " pour écriture." << std::endl;
        return;
    }

    //En-tête du VTK
    fichier << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    fichier << "  <UnstructuredGrid>\n";
    fichier << "    <Piece NumberOfPoints=\"" << nb << "\" NumberOfCells=\"0\">\n";

    //positions des particules
    fichier << "      <Points>\n";
    fichier << "        <DataArray name=\"Position\" type=\"Float32\" NumberOfComponents=\"3\" format=\"binary\">\n";
    for (const auto& cellule : grille) {
        for (const auto& particule : cellule.particules) {
            Vecteur pos = particule.getPosition();
            fichier << std::fixed << std::setprecision(6) << pos[0] << " " << pos[1] << " " << pos[2] << " ";
        }
    }
    fichier << "\n";
    fichier << "        </DataArray>\n";
    fichier << "      </Points>\n";

    //PointData (vitesses et masses des particules)
    fichier << "      <PointData Vectors=\"vector\">\n";
    fichier << "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"binary\">\n";
    for (const auto& cellule : grille) {
        for (const auto& particule : cellule.particules) {
            Vecteur vitesse = particule.getVitesse();
            fichier << fixed << setprecision(6) << vitesse[0] << " " << vitesse[1] << " " << vitesse[2] << " ";
        }
    }
    fichier << "\n";
    fichier << "        </DataArray>\n";

    fichier << "        <DataArray type=\"Float32\" Name=\"Masse\" format=\"binary\">\n";
    for (const auto& cellule : grille) {
        for (const auto& particule : cellule.particules) {
            fichier << fixed << setprecision(6) << particule.getMasse() << " ";
        }
    }
    fichier << "\n";
    fichier << "        </DataArray>\n";
    fichier << "      </PointData>\n";

   
    fichier << "      <Cells>\n";
    fichier << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"binary\">\n";
    fichier << "        </DataArray>\n";
    fichier << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"binary\">\n";
    fichier << "        </DataArray>\n";
    fichier << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"binary\">\n";
    fichier << "        </DataArray>\n";
    fichier << "      </Cells>\n";

    //Fermeture du fichier
    fichier << "    </Piece>\n";
    fichier << "  </UnstructuredGrid>\n";
    fichier << "</VTKFile>\n";

    fichier.close();
    std::cout << "Fichier VTK sauvegardé : " << cheminComplet << std::endl;
}
