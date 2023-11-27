from fltk import cercle, cree_fenetre, mise_a_jour, ferme_fenetre, efface_tout
from fltk import type_ev, donne_ev, attend_ev, texte
from random import randrange, uniform
from math import sqrt, sin, cos, pi
from time import time, sleep
from argparse import ArgumentParser
from ctypes import c_int, c_double, cdll, Structure

WINX = 1500
WINY = 1000
BALLSIZEMIN = 8
BALLSIZEMAX = 20
SPEEDMIN = 60.0
SPEEDMAX = 400.0
FPS = 60.0

'''
typedef struct _ball{
    double size;
    double mass;
    double pos[2];
    double v[2];
} Ball;
'''

class BALL(Structure):
    _fields_ = [('size', c_double),
                ('mass', c_double),
                ('pos', c_double * 2),
                ('v', c_double * 2)]

class Ball:
    def __init__(self):
        self.size = randrange(BALLSIZEMIN, BALLSIZEMAX)
        self.mass = self.size * self.size
        self.pos = (randrange(self.size, WINX - self.size),
                    randrange(self.size, WINY - self.size))
        v = uniform(SPEEDMIN, SPEEDMAX)
        arg = uniform(0, 2 * pi)
        self.v = (v * cos(arg), v * sin(arg))

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

if __name__ == '__main__':
    # Initialization
    parser = ArgumentParser(description='Collision simulator.')
    parser.add_argument('count', metavar='N', type=int, help='the number of balls')
    parser.add_argument('-d', '--diagnostic', action='store_true', 
                        help='whether diagnostics are printed')
    parser.add_argument('-l', '--libpath', nargs='?', type=str,
                        const='libcollision.so', default='libcollision.so',
                        help='relative path of the collision library')
    args = parser.parse_args()    
    cree_fenetre(WINX, WINY)
    # initalize c data
    balls = []
    for _ in range(args.count):
        add_ball(balls)
    cballs = (BALL * args.count)()
    for i in range(args.count):
        cballs[i].size = balls[i].size
        cballs[i].mass = balls[i].mass
        for j in range(2):
            cballs[i].pos[j] = balls[i].pos[j]
            cballs[i].v[j] = balls[i].v[j]
    winsize = (c_double * 2)(WINX, WINY)
    ccnt = c_int(args.count)
    cfps = c_double(1 / FPS)
    # load library
    cmod = cdll.LoadLibrary("./" + args.libpath)
    cmod.update_solid.restype = c_int
    cmod.kinetic_energy.restype = c_double
    cmod.total_momentum.restype = c_double
    # Main iteration
    t1 = time()
    sliding_col = []
    avgcol = 0
    while type_ev(donne_ev()) != 'Quitte':
        efface_tout()
        ret = int(cmod.update_solid(cballs, ccnt, cfps, winsize))
        if ret < 0:
            print("Impossible error!")
            break
        draw_balls(cballs)
        sliding_col.append(ret)
        if len(sliding_col) == 10:
            avgcol = sum(sliding_col) / 10
            sliding_col.clear()
        texte(0, 0, '{:.2f} collisions, avg. in 10 frames'.format(avgcol), taille = 10)
        if args.diagnostic:
            momentum = -1
            energy = float(cmod.kinetic_energy(cballs, ccnt))
            momentum = float(cmod.total_momentum(cballs, ccnt))
            texte(0, 30, 'KE: {:f}\nTM: {:f}'.format(energy, momentum), taille=10)
        tsleep = 1 / FPS - time() + t1
        if tsleep > 0:
            sleep(tsleep)
        else:
            texte(0, 15, 'Late for {:f} ms'.format(-tsleep*1000), taille=10)
        t1 = time()
        mise_a_jour()
