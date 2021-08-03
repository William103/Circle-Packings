from icosanode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs

def functions(n, z):
    if n == 1:
        return -1*np.conj(z)
    duals = [[1/2*(-1 + math.sqrt(5)) - 1j, 1/2*(-1 + math.sqrt(5))], [], 
    [1/22*(1 + 3*math.sqrt(5)) + 1/11*(3 - 2*math.sqrt(5))*1j, 2/(13 + 5*math.sqrt(5))], [1/62*(15 + 13*math.sqrt(5)) + 1/31*(6 - math.sqrt(5))*1j, 
  2/(11 + 7*math.sqrt(5))], [1/22*(1 + 3*math.sqrt(5)) + 
   1/11*(-3 + 2*math.sqrt(5))*1j, 2/(
  13 + 5*math.sqrt(5))], [1/22*(15 - math.sqrt(5)) + 1/11*(-4 + math.sqrt(5))*1j, 2/(
  17 + 7*math.sqrt(5))], [1/22*(15 - math.sqrt(5)) + 1/11*(4 - math.sqrt(5))*1j, 2/(
  17 + 7*math.sqrt(5))], [3/(2*math.sqrt(5)), 1/(
  10 + 4*math.sqrt(5))], [1/math.sqrt(5) + (-(2/5) + 1/math.sqrt(5))*1j, -(2/5) + 1/
   math.sqrt(5)], [-(3/38)*(-9 + math.sqrt(5)) + 1/19*(-8 + 3*math.sqrt(5))*1j, 2/(
  23 + 11*math.sqrt(5))], [-(3/38)*(-9 + math.sqrt(5)) + 1/19*(8 - 3*math.sqrt(5))*1j,
   2/(23 + 11*math.sqrt(5))], [1/6*(-1 + 2*math.sqrt(5)), 
  1/6*(-2 + math.sqrt(5))], [1 + (-2 + math.sqrt(5))*1j, -2 + math.sqrt(5)], [1/(
  2 + 2*math.sqrt(5)), 1/(
  2 + 2*math.sqrt(5))], [1/38*(11 + 3*math.sqrt(5)) + 1/19*(-1 - 2*math.sqrt(5))*1j, 
  2/(7 + 5*math.sqrt(5))], [1/math.sqrt(5) + (2/5 - 1/math.sqrt(5))*1j, -(2/5) + 1/
   math.sqrt(5)], [2/(1 + math.sqrt(5)) + 1j, 2/(
  1 + math.sqrt(5))], [1 + (2 - math.sqrt(5))*1j, -2 + math.sqrt(5)], [1/38*(11 + 3*math.sqrt(5)) + 1/19*(1 + 2*math.sqrt(5))*1j, 2/(
  7 + 5*math.sqrt(5))], [1/62*(15 + 13*math.sqrt(5)) + 1/31*(-6 + math.sqrt(5))*1j, 
  2/(11 + 7*math.sqrt(5))]]
    return duals[n][0] + np.conj((duals[n][1])**2 / (z-duals[n][0]))

def derivatives(n, z):
    if n == 1:
        return 1
    duals = [[1/2*(-1 + math.sqrt(5)) - 1j, 1/2*(-1 + math.sqrt(5))], [], 
    [1/22*(1 + 3*math.sqrt(5)) + 1/11*(3 - 2*math.sqrt(5))*1j, 2/(13 + 5*math.sqrt(5))], [1/62*(15 + 13*math.sqrt(5)) + 1/31*(6 - math.sqrt(5))*1j, 
  2/(11 + 7*math.sqrt(5))], [1/22*(1 + 3*math.sqrt(5)) + 
   1/11*(-3 + 2*math.sqrt(5))*1j, 2/(
  13 + 5*math.sqrt(5))], [1/22*(15 - math.sqrt(5)) + 1/11*(-4 + math.sqrt(5))*1j, 2/(
  17 + 7*math.sqrt(5))], [1/22*(15 - math.sqrt(5)) + 1/11*(4 - math.sqrt(5))*1j, 2/(
  17 + 7*math.sqrt(5))], [3/(2*math.sqrt(5)), 1/(
  10 + 4*math.sqrt(5))], [1/math.sqrt(5) + (-(2/5) + 1/math.sqrt(5))*1j, -(2/5) + 1/
   math.sqrt(5)], [-(3/38)*(-9 + math.sqrt(5)) + 1/19*(-8 + 3*math.sqrt(5))*1j, 2/(
  23 + 11*math.sqrt(5))], [-(3/38)*(-9 + math.sqrt(5)) + 1/19*(8 - 3*math.sqrt(5))*1j,
   2/(23 + 11*math.sqrt(5))], [1/6*(-1 + 2*math.sqrt(5)), 
  1/6*(-2 + math.sqrt(5))], [1 + (-2 + math.sqrt(5))*1j, -2 + math.sqrt(5)], [1/(
  2 + 2*math.sqrt(5)), 1/(
  2 + 2*math.sqrt(5))], [1/38*(11 + 3*math.sqrt(5)) + 1/19*(-1 - 2*math.sqrt(5))*1j, 
  2/(7 + 5*math.sqrt(5))], [1/math.sqrt(5) + (2/5 - 1/math.sqrt(5))*1j, -(2/5) + 1/
   math.sqrt(5)], [2/(1 + math.sqrt(5)) + 1j, 2/(
  1 + math.sqrt(5))], [1 + (2 - math.sqrt(5))*1j, -2 + math.sqrt(5)], [1/38*(11 + 3*math.sqrt(5)) + 1/19*(1 + 2*math.sqrt(5))*1j, 2/(
  7 + 5*math.sqrt(5))], [1/62*(15 + 13*math.sqrt(5)) + 1/31*(-6 + math.sqrt(5))*1j, 
  2/(11 + 7*math.sqrt(5))]]
    return abs(-1*((z-duals[n][0])**2) / ((duals[n][1])**2))

def samplePoint(word):
    points = [2/5*(-1 + math.sqrt(5)) - (4j)/5, -(1/
   2), 1/22*(1 + 3*math.sqrt(5)) + (1/100 + 1/11*(3 - 2*math.sqrt(5)))*1j, 1/
    62*(15 + 13*math.sqrt(5)) + 
   1/31*(6 - math.sqrt(5))*1j, 1/
    22*(1 + 3*math.sqrt(5)) + (-(1/100) + 1/11*(-3 + 2*math.sqrt(5)))*1j, 1/
    22*(15 - math.sqrt(5)) + 
   1/11*(-4 + math.sqrt(5))*1j, 1/22*(15 - math.sqrt(5)) + 
   1/11*(4 - math.sqrt(5))*1j, 3/(
  2*math.sqrt(5)), 1/math.sqrt(5) + (-(2/5) + 1/math.sqrt(5))*1j, -(3/38)*(-9 + math.sqrt(5)) + 
   1/19*(-8 + 3*math.sqrt(5))*1j, -(3/38)*(-9 + math.sqrt(5)) + 
   1/19*(8 - 3*math.sqrt(5))*1j, 1/6*(-1 + 2*math.sqrt(5)), 9/10 + 
   9/10*(-2 + math.sqrt(5))*1j, 1/(
  2 + 2*math.sqrt(5)), 1/38*(11 + 3*math.sqrt(5)) + 
   1/19*(-1 - 2*math.sqrt(5))*1j, 1/math.sqrt(5) + (2/5 - 1/math.sqrt(5))*1j, 8/(
   5*(1 + math.sqrt(5))) + (4j)/5, 9/10 + 
   9/10*(2 - math.sqrt(5))*1j, 1/38*(11 + 3*math.sqrt(5)) + 
   1/19*(1 + 2*math.sqrt(5))*1j, 1/62*(15 + 13*math.sqrt(5)) + 
   1/31*(-6 + math.sqrt(5))*1j]
    p = points[word[-1]]
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    tt = time.time()
    generators = [np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0], [5 + 2*math.sqrt(5), 5 + 2*math.sqrt(5), 0, 0, 0, 0, 0, -1, 
   2*(3 + math.sqrt(5)), 0, 0, 0], [1/2*(7 + math.sqrt(5)), 1/2*(3 + math.sqrt(5)), 
   0, 0, 0, 0, 0, 1/2*(-3 + math.sqrt(5)), 1/2*(7 + math.sqrt(5)), 0, 0, 
   0], [1/2*(5 + 3*math.sqrt(5)), 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, 0, 
   1/2*(1 - math.sqrt(5)), 1/2*(5 + 3*math.sqrt(5)), 0, 0, 
   0], [1/2*(7 + math.sqrt(5)), 2 + math.sqrt(5), 0, 0, 0, 0, 0, 
   1/2*(-3 + math.sqrt(5)), 3, 0, 0, 0], [1/2*(5 + 3*math.sqrt(5)), 3 + math.sqrt(5),
    0, 0, 0, 0, 0, 1/2*(1 - math.sqrt(5)), 3 + 2*math.sqrt(5), 0, 0, 
   0], [2*(3 + math.sqrt(5)), 2*(2 + math.sqrt(5)), 0, 0, 0, 0, 0, -1, 
   2*(3 + math.sqrt(5)), 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
   0], [3 + 2*math.sqrt(5), 3 + math.sqrt(5), 0, 0, 0, 0, 0, 1/2*(1 - math.sqrt(5)), 
   1/2*(5 + 3*math.sqrt(5)), 0, 0, 0], [3, 2 + math.sqrt(5), 0, 0, 0, 0, 0, 
   1/2*(-3 + math.sqrt(5)), 1/2*(7 + math.sqrt(5)), 0, 0, 0], [2*(3 + math.sqrt(5)), 
   5 + 2*math.sqrt(5), 0, 0, 0, 0, 0, -1, 5 + 2*math.sqrt(5), 0, 0, 0]]), np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0], [2*(2 + math.sqrt(5)), 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 
  5 + 2*math.sqrt(5), 0, 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 
  0], [1/2*(7 + 3*math.sqrt(5)), 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 3 + math.sqrt(5),
   0, 0, 0, 0, -1, 0], [1/2*(3 + math.sqrt(5)), 3 + math.sqrt(5), 0, 0, 0, 
  2 + math.sqrt(5), 0, 0, 0, 0, 1/2*(1 - math.sqrt(5)), 0], [0, 0, 0, 0, 0, 1, 
  0, 0, 0, 0, 0, 0], [5 + 2*math.sqrt(5), 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 
  2*(2 + math.sqrt(5)), 0, 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 0], [5 + 2*math.sqrt(5),
   1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 5 + 2*math.sqrt(5), 0, 0, 0, 0, 
  1/2*(-1 - math.sqrt(5)), 0], [2 + math.sqrt(5), 3 + math.sqrt(5), 0, 0, 0, 
  1/2*(3 + math.sqrt(5)), 0, 0, 0, 0, 1/2*(1 - math.sqrt(5)), 0], [2 + math.sqrt(5), 
  1/2*(5 + math.sqrt(5)), 0, 0, 0, 2 + math.sqrt(5), 0, 0, 0, 0, 
  1/2*(1 - math.sqrt(5)), 0], [3 + math.sqrt(5), 2*(2 + math.sqrt(5)), 0, 0, 0, 
  3 + math.sqrt(5), 0, 0, 0, 0, -1, 0], [3 + math.sqrt(5), 1/2*(7 + 3*math.sqrt(5)), 
  0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, -1, 0]]), np.array([[0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 1/2*(5 + 3*math.sqrt(5)), 0, 0, 
   1/2*(1 - math.sqrt(5)), 0, 0, 1/2*(5 + 3*math.sqrt(5)), 0], [0, 1, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0], [0, 1/2*(3 + math.sqrt(5)), 0, 0, 
   1/2*(7 + math.sqrt(5)), 0, 0, 1/2*(-3 + math.sqrt(5)), 0, 0, 
   1/2*(7 + math.sqrt(5)), 0], [0, 5 + 2*math.sqrt(5), 0, 0, 5 + 2*math.sqrt(5), 0, 
   0, -1, 0, 0, 2*(3 + math.sqrt(5)), 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0], [0, 2 + math.sqrt(5), 0, 0, 1/2*(7 + math.sqrt(5)), 0, 0, 
   1/2*(-3 + math.sqrt(5)), 0, 0, 3, 0], [0, 3 + math.sqrt(5), 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 0, 1/2*(1 - math.sqrt(5)), 0, 0, 3 + 2*math.sqrt(5), 
   0], [0, 2*(2 + math.sqrt(5)), 0, 0, 2*(3 + math.sqrt(5)), 0, 0, -1, 0, 0, 
   2*(3 + math.sqrt(5)), 0], [0, 2 + math.sqrt(5), 0, 0, 3, 0, 0, 
   1/2*(-3 + math.sqrt(5)), 0, 0, 1/2*(7 + math.sqrt(5)), 0], [0, 5 + 2*math.sqrt(5),
    0, 0, 2*(3 + math.sqrt(5)), 0, 0, -1, 0, 0, 5 + 2*math.sqrt(5), 0], [0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 3 + math.sqrt(5), 0, 0, 
   3 + 2*math.sqrt(5), 0, 0, 1/2*(1 - math.sqrt(5)), 0, 0, 1/2*(5 + 3*math.sqrt(5)), 
   0]]), np.array([[0, 0, 0, 2 + math.sqrt(5), 0, 1/2*(1 - math.sqrt(5)), 0, 
   1/2*(3 + math.sqrt(5)), 0, 3 + math.sqrt(5), 0, 0], [0, 0, 0, 5 + 2*math.sqrt(5), 
   0, 1/2*(-1 - math.sqrt(5)), 0, 2*(2 + math.sqrt(5)), 0, 1/2*(13 + 5*math.sqrt(5)),
    0, 0], [0, 0, 0, 3 + math.sqrt(5), 0, -1, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
   0], [0, 0, 0, 2*(2 + math.sqrt(5)), 0, 1/2*(-1 - math.sqrt(5)), 0, 
   5 + 2*math.sqrt(5), 0, 1/2*(13 + 5*math.sqrt(5)), 0, 0], [0, 0, 0, 
   3 + math.sqrt(5), 0, -1, 0, 3 + math.sqrt(5), 0, 2*(2 + math.sqrt(5)), 0, 0], [0, 
   0, 0, 2 + math.sqrt(5), 0, 1/2*(1 - math.sqrt(5)), 0, 2 + math.sqrt(5), 0, 
   1/2*(5 + math.sqrt(5)), 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0,
    0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, -1, 0, 3 + math.sqrt(5), 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 
   0], [0, 0, 0, 5 + 2*math.sqrt(5), 0, 1/2*(-1 - math.sqrt(5)), 0, 
   5 + 2*math.sqrt(5), 0, 1/2*(11 + 5*math.sqrt(5)), 0, 0], [0, 0, 0, 
   1/2*(3 + math.sqrt(5)), 0, 1/2*(1 - math.sqrt(5)), 0, 2 + math.sqrt(5), 0, 
   3 + math.sqrt(5), 0, 0]]), np.array([[0, 0, 1/2*(-1 - math.sqrt(5)), 0, 
   1/2*(7 + 3*math.sqrt(5)), 1/2*(5 + 3*math.sqrt(5)), 0, 0, 0, 0, 0, 
   1/2*(7 + 3*math.sqrt(5))], [0, 0, -1, 0, 3 + math.sqrt(5), 1/2*(5 + math.sqrt(5)),
    0, 0, 0, 0, 0, 1/2*(5 + math.sqrt(5))], [0, 0, -1, 0, 3 + math.sqrt(5), 2, 
   0, 0, 0, 0, 0, 3 + math.sqrt(5)], [0, 0, 1/2*(-3 - math.sqrt(5)), 0, 
   1/2*(11 + 5*math.sqrt(5)), 3/2*(3 + math.sqrt(5)), 0, 0, 0, 0, 0, 
   1/2*(13 + 5*math.sqrt(5))], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1/2*(-3 - math.sqrt(5)), 0, 
   1/2*(13 + 5*math.sqrt(5)), 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, 0, 
   1/2*(13 + 5*math.sqrt(5))], [0, 0, 1/2*(-1 - math.sqrt(5)), 0, 
   1/2*(7 + 3*math.sqrt(5)), 2 + math.sqrt(5), 0, 0, 0, 0, 0, 
   2*(2 + math.sqrt(5))], [0, 0, 1/2*(-3 - math.sqrt(5)), 0, 
   1/2*(13 + 5*math.sqrt(5)), 3/2*(3 + math.sqrt(5)), 0, 0, 0, 0, 0, 
   1/2*(11 + 5*math.sqrt(5))], [0, 0, -1, 0, 1/2*(5 + math.sqrt(5)), 
   1/2*(5 + math.sqrt(5)), 0, 0, 0, 0, 0, 3 + math.sqrt(5)], [0, 0, 
   1/2*(-1 - math.sqrt(5)), 0, 2*(2 + math.sqrt(5)), 2 + math.sqrt(5), 0, 0, 0, 0, 0,
    1/2*(7 + 3*math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]), np.array([[0, 
   0, 0, 0, 0, 0, 1/2*(5 + 3*math.sqrt(5)), 0, 3 + 2*math.sqrt(5), 
   1/2*(1 - math.sqrt(5)), 3 + math.sqrt(5), 0], [0, 0, 0, 0, 0, 0, 3, 0, 
   1/2*(7 + math.sqrt(5)), 1/2*(-3 + math.sqrt(5)), 2 + math.sqrt(5), 0], [0, 0, 0, 
   0, 0, 0, 1/2*(7 + math.sqrt(5)), 0, 3, 1/2*(-3 + math.sqrt(5)), 2 + math.sqrt(5), 
   0], [0, 0, 0, 0, 0, 0, 1/2*(7 + math.sqrt(5)), 0, 1/2*(7 + math.sqrt(5)), 
   1/2*(-3 + math.sqrt(5)), 1/2*(3 + math.sqrt(5)), 0], [0, 0, 0, 0, 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 1/2*(5 + 3*math.sqrt(5)), 1/2*(1 - math.sqrt(5)), 
   1/2*(7 + 3*math.sqrt(5)), 0], [0, 0, 0, 0, 0, 0, 5 + 2*math.sqrt(5), 0, 
   2*(3 + math.sqrt(5)), -1, 5 + 2*math.sqrt(5), 0], [0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0], [0, 0, 0, 0, 0, 0, 3 + 2*math.sqrt(5), 0, 
   1/2*(5 + 3*math.sqrt(5)), 1/2*(1 - math.sqrt(5)), 3 + math.sqrt(5), 0], [0, 0, 0, 
   0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 2*(3 + math.sqrt(5)), 0, 
   2*(3 + math.sqrt(5)), -1, 2*(2 + math.sqrt(5)), 0], [0, 0, 0, 0, 0, 0, 0, 0, 
   0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 2*(3 + math.sqrt(5)), 0, 
   5 + 2*math.sqrt(5), -1, 5 + 2*math.sqrt(5), 0]]), np.array([[0, 0, 0, 0, 0, 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 1/2*(7 + 3*math.sqrt(5)), 1/2*(1 - math.sqrt(5)), 
   1/2*(5 + 3*math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 5 + 2*math.sqrt(5), 0, 
   5 + 2*math.sqrt(5), -1, 2*(3 + math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 
   1/2*(7 + math.sqrt(5)), 0, 1/2*(3 + math.sqrt(5)), 1/2*(-3 + math.sqrt(5)), 
   1/2*(7 + math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 1/2*(7 + math.sqrt(5)), 0, 
   2 + math.sqrt(5), 1/2*(-3 + math.sqrt(5)), 3], [0, 0, 0, 0, 0, 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 3 + math.sqrt(5), 1/2*(1 - math.sqrt(5)), 
   3 + 2*math.sqrt(5)], [0, 0, 0, 0, 0, 0, 0, 3, 0, 2 + math.sqrt(5), 
   1/2*(-3 + math.sqrt(5)), 1/2*(7 + math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 
   3 + 2*math.sqrt(5), 0, 3 + math.sqrt(5), 1/2*(1 - math.sqrt(5)), 
   1/2*(5 + 3*math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 
   0, 0, 0, 0, 0, 2*(3 + math.sqrt(5)), 0, 5 + 2*math.sqrt(5), -1, 
   5 + 2*math.sqrt(5)], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 
   0, 0, 0, 2*(3 + math.sqrt(5)), 0, 2*(2 + math.sqrt(5)), -1, 
   2*(3 + math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]), np.array([[0, 0, 0,
    2*(2 + math.sqrt(5)), 0, 0, 2 + math.sqrt(5), 1/2*(7 + 3*math.sqrt(5)), 0, 
   1/2*(-1 - math.sqrt(5)), 0, 0], [0, 0, 0, 1/2*(13 + 5*math.sqrt(5)), 0, 0, 
   3/2*(3 + math.sqrt(5)), 1/2*(11 + 5*math.sqrt(5)), 0, 1/2*(-3 - math.sqrt(5)), 0, 
   0], [0, 0, 0, 1/2*(5 + math.sqrt(5)), 0, 0, 1/2*(5 + math.sqrt(5)), 
   3 + math.sqrt(5), 0, -1, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
   0], [0, 0, 0, 1/2*(11 + 5*math.sqrt(5)), 0, 0, 3/2*(3 + math.sqrt(5)), 
   1/2*(13 + 5*math.sqrt(5)), 0, 1/2*(-3 - math.sqrt(5)), 0, 0], [0, 0, 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 1/2*(7 + 3*math.sqrt(5)), 
   1/2*(13 + 5*math.sqrt(5)), 0, 1/2*(-3 - math.sqrt(5)), 0, 0], [0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0,
    0, 3 + math.sqrt(5), 0, 0, 1/2*(5 + math.sqrt(5)), 1/2*(5 + math.sqrt(5)), 0, -1,
    0, 0], [0, 0, 0, 3 + math.sqrt(5), 0, 0, 2, 3 + math.sqrt(5), 0, -1, 0, 
   0], [0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 1/2*(5 + 3*math.sqrt(5)), 
   1/2*(7 + 3*math.sqrt(5)), 0, 1/2*(-1 - math.sqrt(5)), 0, 0], [0, 0, 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 2 + math.sqrt(5), 2*(2 + math.sqrt(5)), 0, 
   1/2*(-1 - math.sqrt(5)), 0, 0]]), np.array([[0, 0, 1/2*(11 + 5*math.sqrt(5)), 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 0, 0, 1/2*(-3 - math.sqrt(5)), 
   3/2*(3 + math.sqrt(5))], [0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 2*(2 + math.sqrt(5)),
    0, 0, 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 2 + math.sqrt(5)], [0, 0, 1, 0, 0, 
   0, 0, 0, 0, 0, 0, 0], [0, 0, 1/2*(13 + 5*math.sqrt(5)), 0, 
   1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 0, 0, 1/2*(-3 - math.sqrt(5)), 
   3/2*(3 + math.sqrt(5))], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 
   1/2*(5 + math.sqrt(5)), 0, 3 + math.sqrt(5), 0, 0, 0, 0, 0, -1, 
   1/2*(5 + math.sqrt(5))], [0, 0, 2*(2 + math.sqrt(5)), 0, 1/2*(7 + 3*math.sqrt(5)),
    0, 0, 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 2 + math.sqrt(5)], [0, 0, 
   3 + math.sqrt(5), 0, 1/2*(5 + math.sqrt(5)), 0, 0, 0, 0, 0, -1, 
   1/2*(5 + math.sqrt(5))], [0, 0, 1/2*(13 + 5*math.sqrt(5)), 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 0, 0, 1/2*(-3 - math.sqrt(5)), 
   1/2*(7 + 3*math.sqrt(5))], [0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 
   1/2*(5 + 3*math.sqrt(5))], [0, 0, 3 + math.sqrt(5), 0, 3 + math.sqrt(5), 0, 0, 0, 
   0, 0, -1, 2], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]), np.array([[0, 
   1/2*(-1 - math.sqrt(5)), 2*(2 + math.sqrt(5)), 0, 0, 0, 5 + 2*math.sqrt(5), 0, 0, 
   0, 1/2*(13 + 5*math.sqrt(5)), 0], [0, -1, 3 + math.sqrt(5), 0, 0, 0, 
   3 + math.sqrt(5), 0, 0, 0, 2*(2 + math.sqrt(5)), 0], [0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0], [0, -1, 3 + math.sqrt(5), 0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0,
    0, 0, 1/2*(7 + 3*math.sqrt(5)), 0], [0, 1/2*(1 - math.sqrt(5)), 2 + math.sqrt(5),
    0, 0, 0, 1/2*(3 + math.sqrt(5)), 0, 0, 0, 3 + math.sqrt(5), 0], [0, 
   1/2*(-1 - math.sqrt(5)), 5 + 2*math.sqrt(5), 0, 0, 0, 2*(2 + math.sqrt(5)), 0, 0, 
   0, 1/2*(13 + 5*math.sqrt(5)), 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
   0], [0, 1/2*(1 - math.sqrt(5)), 2 + math.sqrt(5), 0, 0, 0, 2 + math.sqrt(5), 0, 0,
    0, 1/2*(5 + math.sqrt(5)), 0], [0, 1/2*(1 - math.sqrt(5)), 
   1/2*(3 + math.sqrt(5)), 0, 0, 0, 2 + math.sqrt(5), 0, 0, 0, 3 + math.sqrt(5), 
   0], [0, 1/2*(-1 - math.sqrt(5)), 5 + 2*math.sqrt(5), 0, 0, 0, 5 + 2*math.sqrt(5), 
   0, 0, 0, 1/2*(11 + 5*math.sqrt(5)), 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   1, 0], [0, -1, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 3 + math.sqrt(5), 0, 0, 0, 
   1/2*(7 + 3*math.sqrt(5)), 0]]), np.array([[0, 0, 1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 
   1/2*(-3 - math.sqrt(5)), 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 
   3/2*(3 + math.sqrt(5))], [0, 0, 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 
   1/2*(-3 - math.sqrt(5)), 1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 
   3/2*(3 + math.sqrt(5))], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 2*(2 + math.sqrt(5)), 
   0, 0, 0, 2 + math.sqrt(5)], [0, 0, 3 + math.sqrt(5), 0, 0, 0, -1, 
   1/2*(5 + math.sqrt(5)), 0, 0, 0, 1/2*(5 + math.sqrt(5))], [0, 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 1/2*(5 + 3*math.sqrt(5))], [0, 0, 
   3 + math.sqrt(5), 0, 0, 0, -1, 3 + math.sqrt(5), 0, 0, 0, 2], [0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0], [0, 0, 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 
   1/2*(-3 - math.sqrt(5)), 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 
   1/2*(7 + 3*math.sqrt(5))], [0, 0, 1/2*(5 + math.sqrt(5)), 0, 0, 0, -1, 
   3 + math.sqrt(5), 0, 0, 0, 1/2*(5 + math.sqrt(5))], [0, 0, 2*(2 + math.sqrt(5)), 
   0, 0, 0, 1/2*(-1 - math.sqrt(5)), 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 
   2 + math.sqrt(5)], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]), np.array([[-1, 0, 
   2*(2 + math.sqrt(5)), 0, 0, 0, 2*(3 + math.sqrt(5)), 2*(3 + math.sqrt(5)), 0, 0, 
   0, 0], [-1, 0, 5 + 2*math.sqrt(5), 0, 0, 0, 2*(3 + math.sqrt(5)), 
   5 + 2*math.sqrt(5), 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 
   0], [1/2*(-3 + math.sqrt(5)), 0, 1/2*(3 + math.sqrt(5)), 0, 0, 0, 
   1/2*(7 + math.sqrt(5)), 1/2*(7 + math.sqrt(5)), 0, 0, 0, 
   0], [1/2*(1 - math.sqrt(5)), 0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 1/2*(5 + 3*math.sqrt(5)), 0, 0, 0, 0], [-1, 0, 
   5 + 2*math.sqrt(5), 0, 0, 0, 5 + 2*math.sqrt(5), 2*(3 + math.sqrt(5)), 0, 0, 0, 
   0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 
   0, 0, 0, 0], [1/2*(1 - math.sqrt(5)), 0, 3 + math.sqrt(5), 0, 0, 0, 
   3 + 2*math.sqrt(5), 1/2*(5 + 3*math.sqrt(5)), 0, 0, 0, 
   0], [1/2*(1 - math.sqrt(5)), 0, 3 + math.sqrt(5), 0, 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 3 + 2*math.sqrt(5), 0, 0, 0, 
   0], [1/2*(-3 + math.sqrt(5)), 0, 2 + math.sqrt(5), 0, 0, 0, 
   1/2*(7 + math.sqrt(5)), 3, 0, 0, 0, 0], [1/2*(-3 + math.sqrt(5)), 0, 
   2 + math.sqrt(5), 0, 0, 0, 3, 1/2*(7 + math.sqrt(5)), 0, 0, 0, 0]]), np.array([[1, 0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1/2*(7 + 3*math.sqrt(5)), 0, 0, 
   3 + math.sqrt(5), 0, 0, 0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 
   0, -1], [2*(2 + math.sqrt(5)), 0, 0, 5 + 2*math.sqrt(5), 0, 0, 0, 0, 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 1/2*(-1 - math.sqrt(5))], [0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0], [5 + 2*math.sqrt(5), 0, 0, 2*(2 + math.sqrt(5)), 0, 0, 0, 0,
    0, 1/2*(13 + 5*math.sqrt(5)), 0, 1/2*(-1 - math.sqrt(5))], [2 + math.sqrt(5), 0, 
   0, 1/2*(3 + math.sqrt(5)), 0, 0, 0, 0, 0, 3 + math.sqrt(5), 0, 
   1/2*(1 - math.sqrt(5))], [3 + math.sqrt(5), 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 
   0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, -1], [1/2*(3 + math.sqrt(5)), 0, 0, 
   2 + math.sqrt(5), 0, 0, 0, 0, 0, 3 + math.sqrt(5), 0, 
   1/2*(1 - math.sqrt(5))], [2 + math.sqrt(5), 0, 0, 2 + math.sqrt(5), 0, 0, 0, 0, 0,
    1/2*(5 + math.sqrt(5)), 0, 1/2*(1 - math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0], [5 + 2*math.sqrt(5), 0, 0, 5 + 2*math.sqrt(5), 0, 0, 0, 0, 0, 
   1/2*(11 + 5*math.sqrt(5)), 0, 1/2*(-1 - math.sqrt(5))], [3 + math.sqrt(5), 0, 0, 
   3 + math.sqrt(5), 0, 0, 0, 0, 0, 2*(2 + math.sqrt(5)), 0, -1]]), np.array([[0, 
   2 + math.sqrt(5), 0, 0, 1/2*(3 + math.sqrt(5)), 3 + math.sqrt(5), 0, 0, 0, 
   1/2*(1 - math.sqrt(5)), 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0,
    3 + math.sqrt(5), 0, 0, 1/2*(7 + 3*math.sqrt(5)), 1/2*(7 + 3*math.sqrt(5)), 0, 0,
    0, -1, 0, 0], [0, 5 + 2*math.sqrt(5), 0, 0, 2*(2 + math.sqrt(5)), 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 0, 0], [0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
   0], [0, 5 + 2*math.sqrt(5), 0, 0, 5 + 2*math.sqrt(5), 1/2*(11 + 5*math.sqrt(5)), 
   0, 0, 0, 1/2*(-1 - math.sqrt(5)), 0, 0], [0, 2*(2 + math.sqrt(5)), 0, 0, 
   5 + 2*math.sqrt(5), 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 1/2*(-1 - math.sqrt(5)), 
   0, 0], [0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 3 + math.sqrt(5), 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, -1, 0, 0], [0, 3 + math.sqrt(5), 0, 0, 
   3 + math.sqrt(5), 2*(2 + math.sqrt(5)), 0, 0, 0, -1, 0, 0], [0, 2 + math.sqrt(5), 
   0, 0, 2 + math.sqrt(5), 1/2*(5 + math.sqrt(5)), 0, 0, 0, 1/2*(1 - math.sqrt(5)), 
   0, 0], [0, 1/2*(3 + math.sqrt(5)), 0, 0, 2 + math.sqrt(5), 3 + math.sqrt(5), 0, 0,
    0, 1/2*(1 - math.sqrt(5)), 0, 0]]), np.array([[0, 3 + math.sqrt(5), 0, 0, 0, 
   1/2*(1 - math.sqrt(5)), 0, 0, 2 + math.sqrt(5), 0, 1/2*(3 + math.sqrt(5)), 0], [0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 
   0, -1, 0, 0, 3 + math.sqrt(5), 0, 1/2*(7 + 3*math.sqrt(5)), 0], [0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, -1, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 
   3 + math.sqrt(5), 0], [0, 3 + math.sqrt(5), 0, 0, 0, 1/2*(1 - math.sqrt(5)), 0, 0,
    1/2*(3 + math.sqrt(5)), 0, 2 + math.sqrt(5), 0], [0, 2*(2 + math.sqrt(5)), 0, 0, 
   0, -1, 0, 0, 3 + math.sqrt(5), 0, 3 + math.sqrt(5), 0], [0, 
   1/2*(5 + math.sqrt(5)), 0, 0, 0, 1/2*(1 - math.sqrt(5)), 0, 0, 2 + math.sqrt(5), 
   0, 2 + math.sqrt(5), 0], [0, 1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 
   1/2*(-1 - math.sqrt(5)), 0, 0, 5 + 2*math.sqrt(5), 0, 5 + 2*math.sqrt(5), 0], [0, 
   0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 1/2*(13 + 5*math.sqrt(5)), 0, 0, 
   0, 1/2*(-1 - math.sqrt(5)), 0, 0, 5 + 2*math.sqrt(5), 0, 2*(2 + math.sqrt(5)), 
   0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 1/2*(13 + 5*math.sqrt(5)),
    0, 0, 0, 1/2*(-1 - math.sqrt(5)), 0, 0, 2*(2 + math.sqrt(5)), 0, 
   5 + 2*math.sqrt(5), 0]]), np.array([[0, 1/2*(-3 - math.sqrt(5)), 1/2*(7 + 3*math.sqrt(5)), 
   0, 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 0, 0, 1/2*(13 + 5*math.sqrt(5)), 
   0], [0, -1, 2, 0, 3 + math.sqrt(5), 0, 0, 0, 0, 0, 3 + math.sqrt(5), 0], [0, 
   0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1/2*(-3 - math.sqrt(5)), 
   3/2*(3 + math.sqrt(5)), 0, 1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 0, 0, 
   1/2*(13 + 5*math.sqrt(5)), 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0,
    1/2*(-1 - math.sqrt(5)), 2 + math.sqrt(5), 0, 2*(2 + math.sqrt(5)), 0, 0, 0, 0, 
   0, 1/2*(7 + 3*math.sqrt(5)), 0], [0, -1, 1/2*(5 + math.sqrt(5)), 0, 
   1/2*(5 + math.sqrt(5)), 0, 0, 0, 0, 0, 3 + math.sqrt(5), 0], [0, 
   1/2*(-1 - math.sqrt(5)), 1/2*(5 + 3*math.sqrt(5)), 0, 1/2*(7 + 3*math.sqrt(5)), 0,
    0, 0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0], [0, 1/2*(-1 - math.sqrt(5)), 
   2 + math.sqrt(5), 0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, 0, 
   2*(2 + math.sqrt(5)), 0], [0, 1/2*(-3 - math.sqrt(5)), 3/2*(3 + math.sqrt(5)), 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 0, 0, 1/2*(11 + 5*math.sqrt(5)), 0], [0, 
   0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, -1, 1/2*(5 + math.sqrt(5)), 0, 
   3 + math.sqrt(5), 0, 0, 0, 0, 0, 1/2*(5 + math.sqrt(5)), 0]]), np.array([[1, 0, 0, 0, 
   0, 0, 0, 0, 0, 0, 0, 0], [2 + math.sqrt(5), 0, 1/2*(-3 + math.sqrt(5)), 0, 0,
    1/2*(7 + math.sqrt(5)), 0, 0, 0, 3, 0, 0], [2*(2 + math.sqrt(5)), 0, -1, 0, 
   0, 2*(3 + math.sqrt(5)), 0, 0, 0, 2*(3 + math.sqrt(5)), 0, 0], [2 + math.sqrt(5), 
   0, 1/2*(-3 + math.sqrt(5)), 0, 0, 3, 0, 0, 0, 1/2*(7 + math.sqrt(5)), 0, 
   0], [3 + math.sqrt(5), 0, 1/2*(1 - math.sqrt(5)), 0, 0, 3 + 2*math.sqrt(5), 0, 0, 
   0, 1/2*(5 + 3*math.sqrt(5)), 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
   0], [5 + 2*math.sqrt(5), 0, -1, 0, 0, 5 + 2*math.sqrt(5), 0, 0, 0, 
   2*(3 + math.sqrt(5)), 0, 0], [3 + math.sqrt(5), 0, 1/2*(1 - math.sqrt(5)), 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 0, 0, 3 + 2*math.sqrt(5), 0, 
   0], [1/2*(7 + 3*math.sqrt(5)), 0, 1/2*(1 - math.sqrt(5)), 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 0, 0, 1/2*(5 + 3*math.sqrt(5)), 0, 0], [0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0], [5 + 2*math.sqrt(5), 0, -1, 0, 0, 
   2*(3 + math.sqrt(5)), 0, 0, 0, 5 + 2*math.sqrt(5), 0, 0], [1/2*(3 + math.sqrt(5)),
    0, 1/2*(-3 + math.sqrt(5)), 0, 0, 1/2*(7 + math.sqrt(5)), 0, 0, 0, 
   1/2*(7 + math.sqrt(5)), 0, 0]]), np.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
   0], [3 + math.sqrt(5), 0, 0, 1/2*(5 + math.sqrt(5)), 0, 0, 0, 0, 
   1/2*(5 + math.sqrt(5)), -1, 0, 0], [1/2*(11 + 5*math.sqrt(5)), 0, 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 0, 3/2*(3 + math.sqrt(5)), 
   1/2*(-3 - math.sqrt(5)), 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
   0], [1/2*(13 + 5*math.sqrt(5)), 0, 0, 1/2*(11 + 5*math.sqrt(5)), 0, 0, 0, 0, 
   3/2*(3 + math.sqrt(5)), 1/2*(-3 - math.sqrt(5)), 0, 0], [2*(2 + math.sqrt(5)), 0, 
   0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, 2 + math.sqrt(5), 
   1/2*(-1 - math.sqrt(5)), 0, 0], [1/2*(5 + math.sqrt(5)), 0, 0, 3 + math.sqrt(5), 
   0, 0, 0, 0, 1/2*(5 + math.sqrt(5)), -1, 0, 0], [1/2*(7 + 3*math.sqrt(5)), 0, 
   0, 2*(2 + math.sqrt(5)), 0, 0, 0, 0, 2 + math.sqrt(5), 1/2*(-1 - math.sqrt(5)), 0,
    0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [3 + math.sqrt(5), 0, 0, 
   3 + math.sqrt(5), 0, 0, 0, 0, 2, -1, 0, 0], [1/2*(7 + 3*math.sqrt(5)), 0, 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0, 0, 1/2*(5 + 3*math.sqrt(5)), 
   1/2*(-1 - math.sqrt(5)), 0, 0], [1/2*(13 + 5*math.sqrt(5)), 0, 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 
   1/2*(-3 - math.sqrt(5)), 0, 0]]), np.array([[0, 0, 0, 0, 0, 1/2*(7 + math.sqrt(5)), 0, 
   0, 0, 2 + math.sqrt(5), 1/2*(-3 + math.sqrt(5)), 3], [0, 0, 0, 0, 0, 
   3 + 2*math.sqrt(5), 0, 0, 0, 3 + math.sqrt(5), 1/2*(1 - math.sqrt(5)), 
   1/2*(5 + 3*math.sqrt(5))], [0, 0, 0, 0, 0, 1/2*(5 + 3*math.sqrt(5)), 0, 0, 0,
    3 + math.sqrt(5), 1/2*(1 - math.sqrt(5)), 3 + 2*math.sqrt(5)], [0, 0, 0, 0, 0, 
   1/2*(5 + 3*math.sqrt(5)), 0, 0, 0, 1/2*(7 + 3*math.sqrt(5)), 
   1/2*(1 - math.sqrt(5)), 1/2*(5 + 3*math.sqrt(5))], [0, 0, 0, 0, 0, 
   1/2*(7 + math.sqrt(5)), 0, 0, 0, 1/2*(3 + math.sqrt(5)), 1/2*(-3 + math.sqrt(5)), 
   1/2*(7 + math.sqrt(5))], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0,
    0, 0, 5 + 2*math.sqrt(5), 0, 0, 0, 5 + 2*math.sqrt(5), -1, 
   2*(3 + math.sqrt(5))], [0, 0, 0, 0, 0, 3, 0, 0, 0, 2 + math.sqrt(5), 
   1/2*(-3 + math.sqrt(5)), 1/2*(7 + math.sqrt(5))], [0, 0, 0, 0, 0, 
   2*(3 + math.sqrt(5)), 0, 0, 0, 5 + 2*math.sqrt(5), -1, 5 + 2*math.sqrt(5)], [0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 2*(3 + math.sqrt(5)), 0,
    0, 0, 2*(2 + math.sqrt(5)), -1, 2*(3 + math.sqrt(5))], [0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1]]), np.array([[0, 1/2*(1 - math.sqrt(5)), 0, 2 + math.sqrt(5), 0, 0, 
   1/2*(3 + math.sqrt(5)), 0, 3 + math.sqrt(5), 0, 0, 0], [0, -1, 0, 
   3 + math.sqrt(5), 0, 0, 3 + math.sqrt(5), 0, 2*(2 + math.sqrt(5)), 0, 0, 
   0], [0, -1, 0, 3 + math.sqrt(5), 0, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 
   0], [0, 1/2*(-1 - math.sqrt(5)), 0, 2*(2 + math.sqrt(5)), 0, 0, 
   5 + 2*math.sqrt(5), 0, 1/2*(13 + 5*math.sqrt(5)), 0, 0, 0], [0, 
   1/2*(-1 - math.sqrt(5)), 0, 5 + 2*math.sqrt(5), 0, 0, 2*(2 + math.sqrt(5)), 0, 
   1/2*(13 + 5*math.sqrt(5)), 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 
   0], [0, 1/2*(1 - math.sqrt(5)), 0, 2 + math.sqrt(5), 0, 0, 2 + math.sqrt(5), 0, 
   1/2*(5 + math.sqrt(5)), 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 
   0], [0, -1, 0, 1/2*(7 + 3*math.sqrt(5)), 0, 0, 3 + math.sqrt(5), 0, 
   1/2*(7 + 3*math.sqrt(5)), 0, 0, 0], [0, 1/2*(1 - math.sqrt(5)), 0, 
   1/2*(3 + math.sqrt(5)), 0, 0, 2 + math.sqrt(5), 0, 3 + math.sqrt(5), 0, 0, 0], [0,
    1/2*(-1 - math.sqrt(5)), 0, 5 + 2*math.sqrt(5), 0, 0, 5 + 2*math.sqrt(5), 0, 
   1/2*(11 + 5*math.sqrt(5)), 0, 0, 0]])]

    root = Node([-1, 2, 10 + 3*math.sqrt(5), 4 + math.sqrt(5), 1/(1 - 2/math.sqrt(5)), 2, 7 + 
 3*math.sqrt(5), 7 + 
 3*math.sqrt(5), 1/2*(7 + math.sqrt(5)), 1/2*(7 + math.sqrt(5)), 1/2*(11 + 
   5*math.sqrt(5)), 1/2*(11 + 5*math.sqrt(5))], [], words, False)

    current_leaves = [root]
    nodes = 1
    while True:
        new_leaves = []
        for leaf in current_leaves:
            next_gen = leaf.next_generation(words, dc, generators)
            new_leaves += next_gen
            nodes += len(next_gen) if len(next_gen) > 1 else 0
        if current_leaves == new_leaves:
            break
        else:
            current_leaves = new_leaves
    for i,leaf in enumerate(current_leaves):
        words[str(leaf.word)] = i
    print(len(current_leaves), "partitions")
    print(nodes,"nodes")
    print("tree construction (s): %f\ntree construction (m): %f" % (time.time()-tt, (time.time()-tt)/60))
    return current_leaves

def constructMatrix(words, dc):
    leave = generateTree(words, dc)
    row = []
    col = []
    data = []
    for i,leaf in enumerate(leave):
        thing = words[str(leaf.word[1:])]
        if isinstance(thing,int):
            row.append(i)
            col.append(thing)
            data.append(sampleValue(leaf.word))
        else:
            sample = sampleValue(leaf.word)
            for wor in thing.leaves():
                row.append(i)
                col.append(words[str(wor.word)])
                data.append(sample)
    return csr_matrix((data,(row,col)),shape=(len(leave),len(leave)))

def secant(x0,y0,x1,y1,z):
    return x0 - (y0-z) * ((x1-x0)/(y1-y0))

def matrixFunction(matrix,l,a):
    #return np.real(eigs(matrix.power(a),k=1)[0][0])
    matrix = matrix.power(a)
    vec = np.ones(l)
    previous_entry = vec[0]
    previous_val = 0
    current = matrix * vec
    current_val = current[0] / previous_entry
    count = 0
    while count < 10000000000 and abs(current_val - previous_val) > 1e-10:
        previous_val = current_val
        previous_entry = current[0]
        current = matrix * current
        #print(current[0],previous_entry)
        current_val = current[0] / previous_entry
        count += 1
    print("power method:", count)
    return current_val

def secantMethod(matrix,l,z,x1,x2,e,its):
    k1 = x1
    k2 = x2
    y1 = matrixFunction(matrix,l,k1)
    y2 = matrixFunction(matrix,l,k2)
    #y1 = testFunction(k1)
    #y2 = testFunction(k2)
    count = 1
    print(count,k1,y1)
    while abs(y1-z)>e and count<its:
        k3 = secant(k1,y1,k2,y2,z)
        k1 = k2
        y1 = y2
        k2 = k3
        y2 = matrixFunction(matrix,l,k2)
        #y2 = testFunction(k2)
        count += 1
        print(count,k1,y1)

def main():
    start = time.time()

    words = {}
    m = 1500
    print("maximum:",m)
    matrix = constructMatrix(words, m)
    #print(matrix)
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*19)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()