#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "collision.h"

/*
* ball->v[0] = x;
* ball->v[1] = y;
*/

static void mulvec(double ball[2], double scalar, double res[2]) {
    res[0] = ball[0] * scalar;
    res[1] = ball[1] * scalar;
}

static void vecadd(double first[2], double second[2], double newVector[2]) {
        newVector[0] = first[0] + second[0];
        newVector[1] = first[1] + second[1];
        return;
}

static void vecdist(double first[2], double second[2], double result[2]) {
    double neg_v2[2] = {-second[0], -second[1]};
    vecadd(first, neg_v2, result);
}

static double vecabs(const double vector[2]) {
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
    double mindist = fabs(proj_dist(relv, relpos));

    if (mindist > cl_r) {
        return -1.0;
    }    
    double tmin = -proj_len(relv, relpos) - sqrt(cl_r * cl_r - mindist * mindist);
    tmin /= vecabs(relv);
    if (tmin < 0) {
        return -1.0;
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
        ytime = (r - y) / vy;
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

static void execute_collide(Ball *b1, Ball *b2) {
    double cvec[2];
    vecdist(b1->pos, b2->pos, cvec);
    double len_cvec = vecabs(cvec);
    double cvec_norm[2] = {cvec[0] / len_cvec, cvec[1] / len_cvec};
    double relv[2];
    vecdist(b1->v, b2->v, relv);
    double relvproj = proj_len(cvec_norm, relv);
    double wr = b2->mass / b1->mass;

    mulvec(cvec_norm, -2 * wr / (1 + wr) * relvproj, cvec_norm);
    vecadd(b1->v, cvec_norm, b1->v);

    mulvec(cvec_norm, 2 / (1 + wr) * relvproj, cvec_norm);
    vecadd(b2->v, cvec_norm, b2->v);
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
    mulvec(ball->v, t, ball->v);
    vecadd(ball->pos, ball->v, ball->pos);
    return;
}

int update_solid(Ball balls[], int count, double elapsed, double winsize[]) {
    double collicnt = 0;
    double tmin = 0;
    double localmin = 0;
    double reftype = 0;
    double ref = 0;
    int i, j, b;
    double machin[2] = {localmin, reftype};
    while(1) {
        tmin = 0;
        double minpair[3] = {0, 0, 0};
        for(i = 0; i < count; ++i) {
            for(j = i; j < count; ++j) {
                if(j != i) {
                    localmin = check_collide(balls[i], balls[j]);
                    reftype = -1;
                } else {
                    check_collide_wall(balls[i], machin);
                }
                if(localmin != 0.0) {
                    if(tmin == 0.0 || tmin < localmin) {
                        tmin = localmin;
                        minpair[0] = i;
                        minpair[1] = j;
                        minpair[2] = reftype;
                    }
                }
            }
        }
        if(tmin == 0.0) {
            fprintf(stderr, "impossible error!\n");
            exit(-1);
        }
        if(tmin > elapsed) {
            break;
        }
        for(b = 0; b < count; ++b) {
            rectilinear(&(balls[b]), tmin);
        }
        double idx1 = minpair[0];
        double idx2 = minpair[1];
        ref = minpair[2];
        if(idx1 != idx2) execute_collide(&(balls[(int) idx1]), &(balls[(int) idx2]));
        else execute_collide_wall(&(balls[(int) idx1]), ref);
        collicnt += 1;
        elapsed -= tmin;
    }
    for(b = 0; b < count; ++b) rectilinear(&balls[b], elapsed);
    return (int) collicnt;
}

/*
def update_solid(balls, count, elapsed):
    # try to spend the elapsed time by executing collisions one at a time
    collicnt = 0
    while True:
        tmin = None
        minpair = None
        # check for the next collision
        for i in range(count):
            for j in range(i, count): # case j = i for wall-checking
                if j != i:
                    localmin = check_collide(balls[i], balls[j])
                    reftype = -1
                else:
                    localmin, reftype = check_collide_wall(balls[i])
                if localmin is not None:
                    if tmin is None or tmin > localmin:
                        tmin = localmin
                        minpair = i, j, reftype
        # there should always be some collision
        if tmin is None:
            print("Impossible error!")
            quit()
        # no collision at all for the remaining time
        if tmin > elapsed:
            break
        for b in balls:
            rectilinear(b, tmin)
        idx1, idx2, ref = minpair
        if idx1 != idx2:
            execute_collide(balls[idx1], balls[idx2])
        else:
            execute_collide_wall(balls[idx1], ref)
        collicnt += 1
        elapsed -= tmin
    for b in balls:
        rectilinear(b, elapsed)
    return collicnt

*/

double kinetic_energy(const Ball balls[], int count) {
    double sum = 0;
    int b;
    for(b = 0; b < count; ++b) {
        sum += balls[b].mass * (pow(vecabs(balls[b].v), 2) / 2);
    }
    return sum;
}

double total_momentum(const Ball balls[], int count) {
    double accu[2] = {0};
    int b;
    for(b = 0; b < count; ++b) {
        mulvec(balls[b].v, balls[b].mass, balls[b].v);
        vecadd(accu, balls[b].v, accu);
    }
    return vecabs(accu);
}
