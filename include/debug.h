#ifndef DEBUG_H
#define DEBUG_H

// Définir DEBUG pour activer les messages de débogage
#define DEBUG

#ifdef DEBUG
#define DEBUG_COUT(x) std::cout << x
#else
#define DEBUG_COUT(x)
#endif

#endif // DEBUG_H