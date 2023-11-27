from fltk import cercle, cree_fenetre, mise_a_jour, ferme_fenetre, efface_tout
from fltk import type_ev, donne_ev, attend_ev, texte
from random import randrange, uniform
from math import sqrt, sin, cos, pi
from time import time, sleep
import argparse

WINX = 1500
WINY = 1000
BALLSIZEMIN = 8
BALLSIZEMAX = 20
SPEEDMIN = 60.0
SPEEDMAX = 400.0
FPS = 60.0

class Ball:
    def __init__(self):
        self.size = randrange(BALLSIZEMIN, BALLSIZEMAX)
        self.mass = self.size * self.size
        self.pos = (randrange(self.size, WINX - self.size),
                    randrange(self.size, WINY - self.size))
        v = uniform(SPEEDMIN, SPEEDMAX)
        arg = uniform(0, 2 * pi)
        self.v = (v * cos(arg), v * sin(arg))

####################
# Computation part #
####################

def mulvec(v, r):
    # multiply by a scalar
    return v[0] * r, v[1] * r

def vecadd(v1, v2):
    # adding two vectors
    return v1[0] + v2[0], v1[1] + v2[1]

def vecdist(v1, v2):
    # vector linking two points expressed as vectors
    return vecadd(v1, mulvec(v2, -1))

def vecabs(v):
    # norm of a vector
    return sqrt(v[0] * v[0] + v[1] * v[1])

def sprod(v1, v2):
    # scarlar product of two vectors
    return v1[0] * v2[0] + v1[1] * v2[1]

def proj_len(direction, v):
    # projection of the vector v on the given direction
    return sprod(direction, v) / vecabs(direction)

def vprod(v1, v2):
    # cross product of two vectors
    return v1[0] * v2[1] - v1[1] * v2[0]

def proj_dist(direction, v):
    # the length of the projection of the vector v on the given direction
    return vprod(direction, v) / vecabs(direction)

def check_collide(b1, b2):
    cl_r = b1.size + b2.size
    # Use relative position and speed for collision detection
    relpos = vecdist(b1.pos, b2.pos)
    relv = vecdist(b1.v, b2.v)
    mindist = abs(proj_dist(relv, relpos))
    if mindist > cl_r:
        return None
    tmin = - proj_len(relv, relpos) - sqrt(cl_r * cl_r - mindist * mindist)
    tmin /= vecabs(relv)
    # If not in the good direction or not in the good interval, then do nothing
    if tmin < 0:
        return None
    return tmin

def check_collide_wall(ball):
    x, y = ball.pos
    vx, vy = ball.v
    r = ball.size
    xtime, ytime = 0, 0
    if vx > 0:
        xtime = (WINX - r - x) / vx
    else:
        xtime = (r - x) / vx
    if vy > 0:
        ytime = (WINY - r - y) / vy
    else:
        ytime = (r - y) / vy
    return (xtime, 0) if xtime < ytime else (ytime, 1)

def execute_collide(b1, b2):
    # compute the new velocity after collision
    cvec = vecdist(b1.pos, b2.pos)
    cvec = mulvec(cvec, 1 / vecabs(cvec))
    relv = vecdist(b1.v, b2.v)
    relvproj = proj_len(cvec, relv)
    wr = b2.mass / b1.mass
    b1.v = vecadd(b1.v, mulvec(cvec, -2 * wr / (1 + wr) * relvproj))
    b2.v = vecadd(b2.v, mulvec(cvec, 2 / (1 + wr) * relvproj))
    return

def execute_collide_wall(ball, ref):
    # compute the new velocity after collision
    vx, vy = ball.v
    # Reflection principle
    if ref == 0:
        vx = -vx
    elif ref == 1:
        vy = -vy
    ball.v = vx, vy
    return

def rectilinear(ball, t):
    # perform rectilinear movement of the ball for a duration t
    ball.pos = vecadd(ball.pos, mulvec(ball.v, t))
    return

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

def kinetic_energy(balls):
    # for debug and diagnostic, computes the kinetic energy of the system
    return sum([b.mass * vecabs(b.v) ** 2 / 2 for b in balls])

def total_momentum(balls):
    # for debug and diagnostic, computes the total momentum of the system
    accu = 0, 0
    for b in balls:
        accu = vecadd(accu, mulvec(b.v, b.mass))
    return vecabs(accu)

###########################
# End of computation part #
###########################

def draw_balls(balls):
    for b in balls:
        cercle(b.pos[0], b.pos[1], b.size, remplissage = "red")
    return

def is_overlapping(b1, b2):
    diffx = b1.pos[0] - b2.pos[0]
    diffy = b1.pos[1] - b2.pos[1]
    dist = sqrt(diffx * diffx + diffy * diffy)
    return dist < b1.size + b2.size

def add_ball(balls):
    newball = Ball()
    # Choose random position until balls are not overlapped
    while any(is_overlapping(newball, b) for b in balls):
        newball = Ball()  
    balls.append(newball)
    return

def print_diagnostics(balls):
    energy = kinetic_energy(balls)
    momentum = total_momentum(balls)
    texte(0, 30, 'KE: {:f}\nTM: {:f}'.format(energy, momentum), taille=10)
    return

if __name__ == '__main__':
    # Initialization
    parser = argparse.ArgumentParser(description='Collision simulator.')
    parser.add_argument('count', metavar='N', type=int, help='the number of balls')
    parser.add_argument('-d', '--diagnostic', action='store_true', 
                        help='whether diagnostics are printed')
    args = parser.parse_args()    
    cree_fenetre(WINX, WINY)
    balls = []
    for _ in range(args.count):
        add_ball(balls)
    # Main iteration
    t1 = time()
    sliding_col = []
    avgcol = 0
    while type_ev(donne_ev()) != 'Quitte':
        efface_tout()
        ret = update_solid(balls, args.count, 1 / FPS)
        draw_balls(balls)
        sliding_col.append(ret)
        if len(sliding_col) == 10:
            avgcol = sum(sliding_col) / 10
            sliding_col.clear()
        texte(0, 0, '{:.2f} collisions, avg. in 10 frames'.format(avgcol), taille = 10)
        if args.diagnostic:
            print_diagnostics(balls)
        tsleep = 1 / FPS - time() + t1
        if tsleep > 0:
            sleep(tsleep)
        else:
            texte(0, 15, 'Late for {:f} ms'.format(-tsleep*1000), taille=10)
        t1 = time()
        mise_a_jour()

