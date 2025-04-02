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

Univers::Univers(int nb, int dimension, vector<double> ld, double rcut)
    : nb(nb), dimension(dimension), ld(ld), rcut(rcut) {
    ncd.resize(dimension);
    int totalCells = 1;
    DEBUG_COUT("Début simulation avec " << grille.size() << " cellules");
    // Calcul des ncd et taille totale de la grille
    for (int i = 0; i < dimension; i++) {
        ncd[i] = static_cast<int>(ld[i] / rcut);
        totalCells *= ncd[i];
    }

    grille.resize(totalCells);

    // Distance minimale entre les particules et les bords de la grille
    double marge = rcut / 2.0;

    // Placement des particules
    for (int i = 0; i < nb; i++) {
        vector<double> pos(dimension);

        //les coordonnées de la particule de manière uniforme avec marges
        for (int d = 0; d < dimension; d++) {
            // Nombre de particules par dimension
            int nbParticulesParDimension = static_cast<int>(pow(nb, 1.0 / dimension));

            //l'indice de la particule
            int indice = (i / static_cast<int>(pow(nbParticulesParDimension, d))) % nbParticulesParDimension;

            // position de la particule
            double espaceDisponible = ld[d] - 2 * marge; 
            double pas = espaceDisponible / (nbParticulesParDimension - 1); // Distance entre les particules
            pos[d] = marge + indice * pas; //position de la particule
        }

        Vecteur position(
            pos[0], 
            (dimension > 1) ? pos[1] : 0.0, 
            (dimension > 2) ? pos[2] : 0.0
        );

        //Création de la particule
        Particule particule(i, position);
        particule.afficher();

        //Déterminer la cellule de la grille où placer la particule
        vector<int> cellCoords(dimension);
        for (int d = 0; d < dimension; d++) {
            cellCoords[d] = static_cast<int>(pos[d] / rcut);
        }

        int cellIndex = getIndice(cellCoords);
        grille[cellIndex].particules.push_back(particule);
    }

    //Initialiser les cellules voisines
    initialiserCellulesVoisines();
}



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

            // Determine the cell index et ajoute la particule à la grille
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



void Univers::forceParticule(double sigma, double epsilon) {
    if (grille.empty()) return;
    // Réinitialisation
    Vecteur forceNulle(0.0, 0.0, 0.0);
    for (auto& cellule : grille) {
        for (auto& p : cellule.particules) {
            p.setForce(forceNulle); // Passage par référence
        }
    }

    //Calcul des forces
    for (size_t i = 0; i < grille.size(); ++i) {
        auto& cellule_i = grille[i];
        
        // Interactions intra-cellule
        for (size_t p1 = 0; p1 < cellule_i.particules.size(); ++p1) {
            for (size_t p2 = p1 + 1; p2 < cellule_i.particules.size(); ++p2) {
                Vecteur f = cellule_i.particules[p1].forceInteraction(cellule_i.particules[p2], sigma, epsilon);
                f += cellule_i.particules[p1].forceInteractionGravitationelle(cellule_i.particules[p2]);

                // Mise à jour symétrique
                Vecteur newForceP1 = cellule_i.particules[p1].getForce() + f;
                Vecteur newForceP2 = cellule_i.particules[p2].getForce() - f;
                cellule_i.particules[p1].setForce(newForceP1); 
                cellule_i.particules[p2].setForce(newForceP2); 
            }
        }

        // Interactions inter-cellules
        for (auto& p : cellule_i.particules) {
            for (int j : cellule_i.voisins) {
                if (j <= static_cast<int>(i)) continue; 

                Vecteur centreVoisin = getCentreCellule(j);
                if ((centreVoisin - p.getPosition()).norm() < rcut) {
                    for (auto& p_voisin : grille[j].particules) {
                        Vecteur delta = p_voisin.getPosition() - p.getPosition();
                        if (delta.norm() < rcut) {
                            Vecteur f = p.forceInteraction(p_voisin, sigma, epsilon);
                            f += p.forceInteractionGravitationelle(p_voisin);

                            //Mettre à jour 
                            Vecteur newForceP = p.getForce() + f;
                            Vecteur newForceVoisin = p_voisin.getForce() - f;
                            p.setForce(newForceP); 
                            p_voisin.setForce(newForceVoisin);
                        }
                    }
                }
            }
        }
    }
}


void Univers::mettreAJourPositionDansLimites(Particule& particule) {
    Vecteur pos = particule.getPosition();

    for (size_t d = 0; d < dimension; d++) {  
        if (ld[d] <= 0) continue;
        pos[d] = fmod(pos[d], ld[d]); 
        if (pos[d] < 0) pos[d] += ld[d]; 
    }

    particule.setPosition(pos);  
};



void Univers::miseAJourCellules() {

    vector<Particule> allParticles;

   for (auto& cellule : grille) {
    allParticles.insert(allParticles.end(), cellule.particules.begin(), cellule.particules.end());
    cellule.particules.clear();
    }

    //Reassign les particules à leurs nouvelles cellules
    for (auto& p : allParticles) {
        vector<int> newCellCoords(dimension);
        bool valid = true;
        
        for (int d = 0; d < dimension; d++) {
            double pos = p.getPosition()[d];
            if (pos < 0 || pos >= ld[d]) {  //validation des limites
                valid = false;
                break;
            }
            newCellCoords[d] = static_cast<int>(pos / rcut);
        }

        if (valid) {
            int cellIndex = getIndice(newCellCoords);
            if (cellIndex >= 0 && cellIndex < grille.size()) {
                grille[cellIndex].ajouterParticule(p);
            }
        }
    }
    initialiserCellulesVoisines();
};



void Univers::etatUnivers(double dt,double sigma, double epsilon,double t_end){
    double t = 0.0;
    int pas=1;
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
    string cheminComplet = string(getenv("HOME")) + "/TP_errabhii_madakcha/" + nomFichier;
    cout << "Chemin du fichier : " << cheminComplet << endl;
    
    ofstream fichier(cheminComplet);
    if (!fichier.is_open()) {
        cerr << "Erreur : Impossible d'ouvrir le fichier " << cheminComplet << endl;
        return;
    }

    // Compter le nombre total de particules
    size_t totalParticles = 0;
    for (const auto& cellule : grille) {
        totalParticles += cellule.particules.size();
    }

    // En-tête VTK
    fichier << std::fixed << std::setprecision(6);
    fichier << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    fichier << "  <UnstructuredGrid>\n";
    fichier << "    <Piece NumberOfPoints=\"" << totalParticles << "\" NumberOfCells=\"0\">\n";

    // Positions
    fichier << "      <Points>\n";
    fichier << "        <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& cellule : grille) {
        for (const auto& particule : cellule.particules) {
            Vecteur pos = particule.getPosition();
            fichier << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
        }
    }
    fichier << "        </DataArray>\n";
    fichier << "      </Points>\n";

    // Données (Vitesse + Masse)
    fichier << "      <PointData Vectors=\"Velocity\">\n";
    
    // Vitesse
    fichier << "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& cellule : grille) {
        for (const auto& particule : cellule.particules) {
            Vecteur vitesse = particule.getVitesse();
            fichier << vitesse[0] << " " << vitesse[1] << " " << vitesse[2] << "\n";
        }
    }
    fichier << "        </DataArray>\n";

    // Masse
    fichier << "        <DataArray type=\"Float32\" Name=\"Masse\" format=\"ascii\">\n";
    for (const auto& cellule : grille) {
        for (const auto& particule : cellule.particules) {
            fichier << particule.getMasse() << "\n";
        }
    }
    fichier << "        </DataArray>\n";
    fichier << "      </PointData>\n";

    // Cellules (vide)
    fichier << "      <Cells>\n";
    fichier << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\"/>\n";
    fichier << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\"/>\n";
    fichier << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\"/>\n";
    fichier << "      </Cells>\n";

    // Fermeture
    fichier << "    </Piece>\n";
    fichier << "  </UnstructuredGrid>\n";
    fichier << "</VTKFile>\n";

    fichier.close();
    cout << "Fichier VTK généré avec succès : " << cheminComplet << endl;
}

