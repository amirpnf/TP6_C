#ifndef __COLLISION_H
#define __COLLISION_H

#include <math.h>

#define WINX 1500
#define WINY 1000
#define BALLSIZEMIN 8
#define BALLSIZEMAX 20
#define SPEEDMIN 60.0
#define SPEEDMAX 400.0
#define FPS 60.0

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

static void mulvec(double ball[2], double scalar, double res[2]);

static void vecadd(double first[2], double second[2], double newVector[2]);

static void vecdist(double first[2], double second[2], double newVector[2]);

static double vecabs(const double vector[2]);

static double sprod(double first[2], double second[2]);

static double proj_len(double direction[2], double vector[2]);

static double vprod(double first[2], double second[2]);

static double proj_dist(double direction[2], double vector[2]);

static double check_collide(Ball first, Ball second);

static void check_collide_wall(Ball ball, double res[2]);

static void execute_collide(Ball* b1, Ball* b2);

static void execute_collide_wall(Ball* ball, int ref);

static void rectilinear(Ball* ball, double t);

int update_solid(Ball balls[], int count, double elapsed, double winsize[]);

double kinetic_energy(const Ball balls[], int count);

double total_momentum(const Ball balls[], int count);
#endif