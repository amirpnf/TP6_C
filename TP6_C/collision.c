#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "collision.h"

/*
* ball->v[0] = x;
* ball->v[1] = y;
*/

#define WINX 1500
#define WINY 1000
#define BALLSIZEMIN 8
#define BALLSIZEMAX 20
#define SPEEDMIN 60.0
#define SPEEDMAX 400.0
#define FPS 60.0

static void mulvec(double ball[2], double scalar) {
    ball[0] *= scalar;
    ball[1] *= scalar;
}

static void vecadd(double first[2], double second[2], double newVector[2]) {
    if(newVector != NULL) {
        newVector[0] = first[0] + second[0];
        newVector[1] = first[1] + second[1];
        return;
    }
    fprintf(stderr, "New Vector is NULL ! \n");
}

static void vecdist(double first[2], double second[2], double newVector[2]) {
    if(newVector != NULL) {
        double tmp[2];
        tmp[0] = first[0];
        tmp[1] = first[1];
        mulvec(tmp, -1);
        vecadd(first, tmp, newVector);
    }
    fprintf(stderr, "New Vector is NULL ! \n");
}

static double vecabs(double vector[2]) {
    return sqrt(vector[0] * vector[0] + vector[1] * vector[1]);
}

static double sprod(double first[2], double second[2]) {
    return first[0] * second[0] + first[1] * second[1];
}

static double proj_len(double direction[2], double vector[2]) {
    return sprod(direction, vector) / vecabs(direction);
}

static double vprod(double first[2], double second[2]) {
    return (first[0] * second[1] - first[1] * second[0]);
}

static double proj_dist(double direction[2], double vector[2]) {
    return vprod(direction, vector) / vecabs(direction);
}

static double check_collide(Ball first, Ball second) {
    double cl_r = first.size + second.size;
    double relpos[2] = {0};
    double relv[2] = {0};
    vecdist(first.pos, second.pos, relpos);
    vecdist(first.v, second.v, relv);
    double mindist = abs(proj_dist(relv, relpos));

    if (mindist > cl_r) {
        return -1;
    }    
    double tmin = -proj_len(relv, relpos) - sqrt(cl_r * cl_r - mindist * mindist);
    tmin /= vecabs(relv);
    // If not in the good direction or not in the good interval, then do nothing
    if (tmin < 0) {
        return -1;
    }    
    return tmin;
}

static void check_collide_wall(Ball ball, double res[2]) {
    double x = ball.pos[0];
    double y = ball.pos[1];
    double vx = ball.v[0];
    double vy = ball.v[1];
    double r = ball.size;
    double xtime = 0, ytime = 0;

    if(vx > 0) {
        xtime = (WINX - r - x) / vx;
    } else {
        xtime = (r - x) / vx;
    }
    if(vy > 0) {
        ytime = (WINY - r - y) / vy;
    } else {
        ytime = (WINY - r - y) / vy;
    }
    if(xtime < ytime) {
        res[0] = xtime;
        res[1] = 0;
        return;
    } else {
        res[0] = ytime;
        res[1] = 1;
    }
}

static void execute_collide(Ball* b1, Ball* b2) {
    double cvec[2] = {0}; 
    vecdist(b1->pos, b2->pos, cvec);
    mulvec(cvec, 1 / vecabs(cvec));
    double relv[2] = {0};
    vecdist(b1->v, b2->v, relv);
    double relvproj = proj_len(cvec, relv);
    double wr = b2->mass / b1->mass;
    mulvec(cvec, -2 * wr / (1 + wr) * relvproj);
    vecadd(b1->v, cvec, b1->v);
    mulvec(cvec, 2 / (1 + wr) * relvproj);
    vecadd(b2->v, cvec, b2->v);
    return;
}

static void execute_collide_wall(Ball* ball, int ref) {
    double vx = ball->v[0];
    double vy = ball->v[1];

    if(ref == 0) {
        vx *= -1;
    } else if(ref == 1) {
        vy *= -1;
    }
    ball->v[0] = vx;
    ball->v[1] = vy;
    return;
}

static void rectilinear(Ball* ball, double t) {
    mulvec(ball->v, t);
    vecadd(ball->pos, ball->v, ball->v);
    return;
}

int main(int argc, char* argv[]) {
    return 0;
}