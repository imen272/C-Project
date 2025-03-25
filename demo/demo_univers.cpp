#include "Univers.h"

int main() {
    //the universe
    Univers univers(1600, 6400, 2, {250.0, 40.0}, 2.5 * 1.0, 1.0);

    //Run the simulation until t = 19.5
    univers.etatUnivers(19.5, 1.0, 5.0,19.5);
    
    return 0;
}
