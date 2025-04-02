#include "Univers.h"

int main() {
    //the universe
    Univers univers(1600, 6400, 2, {250.0, 40.0}, 2.5 * 1.0, 1.0);
    univers.sauvegarderVTK("output/particules_t_0.vtu");
    //Run the simulation until t = 19.5
    univers.etatUnivers(1, 1.0, 5.0,19);
    
    return 0;
}
