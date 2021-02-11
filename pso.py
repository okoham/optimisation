from copy import copy
import numpy as np
import random

# Kochenderfer/Wheeler, Algorithm 9.12
class Particle(object):
    # x: array, v, x_best
    pass 

def particle_swarm_optimisation(f, population, k_max, w=1, c1=1, c2=1):
    n = len(population[0].x)
    x_best = copy(population[0].x_best)
    y_best = np.inf
    for p in population:
        y = f(p.x)
        if y < y_best:
            x_best = P.x
            y_best = y
    for k in range(k_max):
        for p in population:
            r1 = random.random()
            r2 = random.random()
            p.x += p.v 
            p.v = w*p.v + c1*r1 * (p.x_best - p.x) + c2*r2* (x_best - p.x)
            y = f(p.x)
            if y < y_best:
                x_best = p.x
                y_best = y 
            if y < f(p.x_best):
                p.x_best = p.x 
    return population 

