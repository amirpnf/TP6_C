#ifndef __COLLISION_H
#define __COLLISION_H

#include <math.h>

typedef struct _ball{
    double size;
    double mass;
    double pos[2];
    double v[2];
} Ball;

/* 
    Main function, update the system for elapsed time
    returns the number of collisions in the duration, returns -1 when error 
*/
int update_solid(Ball balls[], int count, double elapsed, double winsize[]);

/* Diagnostic information, kinetic energy of the system, should be conserved */
double kinetic_energy(const Ball balls[], int count);

/* Diagnostic information, total momentum of the system */
double total_momentum(const Ball balls[], int count);

#endif