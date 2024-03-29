from octanode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix

def functions(n, z):
    """
    funs = [[0.446462692171690 + 0.184930902191118j, 
  0.184930902191118], [0.184930902191118 + 0.446462692171690j, 
  0.184930902191118], [-0.184930902191118 + 0.446462692171690j, 
  0.184930902191118], [-0.446462692171690 + 0.184930902191118j, 
  0.184930902191118], [-0.446462692171690 - 0.184930902191118j, 
  0.184930902191118], [-0.184930902191118 - 0.446462692171690j, 
  0.184930902191118], [0.184930902191118 - 0.446462692171690j, 
  0.184930902191118], [0.446462692171690 - 0.184930902191118j, 
  0.184930902191118], [1.00000000000000 + 0.414213562373095j, 
  0.414213562373095], [0.414213562373095 + 1.00000000000000j, 
  0.414213562373095], [-0.414213562373095 + 1.00000000000000j, 
  0.414213562373095], [-1.00000000000000 + 0.414213562373095j, 
  0.414213562373095], [-1.00000000000000 - 0.414213562373095j, 
  0.414213562373095], [-0.414213562373095 - 1.00000000000000j, 
  0.414213562373095], [0.414213562373095 - 1.00000000000000j, 
  0.414213562373095], [1.00000000000000 - 0.414213562373095j, 
  0.414213562373095]]
    return funs[n][0] + np.conj((funs[n][1])**2 / (z-funs[n][1]))
    """
    if n == 0:
        return 0.446462692171690 + 0.184930902191118j + np.conj(0.0341994385852209/(-0.446462692171690 - 0.184930902191118j + z)) 
    elif n == 1:
        return 0.184930902191118 + 0.446462692171690j + np.conj(0.0341994385852209/(-0.184930902191118 - 0.446462692171690j + z)) 
    elif n == 2:
        return -0.184930902191118 + 0.446462692171690j + np.conj(0.0341994385852209/(0.184930902191118 - 0.446462692171690j + z)) 
    elif n == 3:
        return -0.446462692171690 + 0.184930902191118j + np.conj(0.0341994385852209/(0.446462692171690 - 0.184930902191118j + z)) 
    elif n == 4:
        return -0.446462692171690 - 0.184930902191118j + np.conj(0.0341994385852209/(0.446462692171690 + 0.184930902191118j + z))
    elif n == 5:
        return -0.184930902191118 - 0.446462692171690j + np.conj(0.0341994385852209/(0.184930902191118 + 0.446462692171690j + z))
    elif n == 6:
        return 0.184930902191118 - 0.446462692171690j + np.conj(0.0341994385852209/(-0.184930902191118 + 0.446462692171690j + z))
    elif n == 7:
        return 0.446462692171690 - 0.184930902191118j + np.conj(0.0341994385852209/(-0.446462692171690 + 0.184930902191118j + z))
    elif n == 8:
        return 1.00000000000000 + 0.414213562373095j + np.conj(0.171572875253810/(-1.00000000000000 - 0.414213562373095j + z))
    elif n == 9:
        return 0.414213562373095 + 1.00000000000000j + np.conj(0.171572875253810/(-0.414213562373095 - 1.00000000000000j + z))
    elif n == 10:
        return -0.414213562373095 + 1.00000000000000j + np.conj(0.171572875253810/(0.414213562373095 - 1.00000000000000j + z))
    elif n == 11:
        return -1.00000000000000 + 0.414213562373095j + np.conj(0.171572875253810/(1.00000000000000 - 0.414213562373095j + z))
    elif n == 12:
        return -1.00000000000000 - 0.414213562373095j + np.conj(0.171572875253810/(1.00000000000000 + 0.414213562373095j + z))
    elif n == 13:
        return  -0.414213562373095 - 1.00000000000000j + np.conj(0.171572875253810/(0.414213562373095 + 1.00000000000000j + z))
    elif n == 14:
        return  0.414213562373095 - 1.00000000000000j + np.conj(0.171572875253810/(-0.414213562373095 + 1.00000000000000j + z))
    elif n == 15:
        return 1.00000000000000 - 0.414213562373095j + np.conj(0.171572875253810/(-1.00000000000000 + 0.414213562373095j + z))

def derivatives(n, z):
    """
    funs = [[0.446462692171690 + 0.184930902191118j, 
  0.184930902191118], [0.184930902191118 + 0.446462692171690j, 
  0.184930902191118], [-0.184930902191118 + 0.446462692171690j, 
  0.184930902191118], [-0.446462692171690 + 0.184930902191118j, 
  0.184930902191118], [-0.446462692171690 - 0.184930902191118j, 
  0.184930902191118], [-0.184930902191118 - 0.446462692171690j, 
  0.184930902191118], [0.184930902191118 - 0.446462692171690j, 
  0.184930902191118], [0.446462692171690 - 0.184930902191118j, 
  0.184930902191118], [1.00000000000000 + 0.414213562373095j, 
  0.414213562373095], [0.414213562373095 + 1.00000000000000j, 
  0.414213562373095], [-0.414213562373095 + 1.00000000000000j, 
  0.414213562373095], [-1.00000000000000 + 0.414213562373095j, 
  0.414213562373095], [-1.00000000000000 - 0.414213562373095j, 
  0.414213562373095], [-0.414213562373095 - 1.00000000000000j, 
  0.414213562373095], [0.414213562373095 - 1.00000000000000j, 
  0.414213562373095], [1.00000000000000 - 0.414213562373095j, 
  0.414213562373095]]
    return abs(-1*((z-funs[n][0])**2) / (funs[n][1]**2))
    """
    if n == 0:
        return abs(-29.2402460791314*(-0.446462692171690 - 0.184930902191118j + z)**2)
    elif n == 1:
        return abs(-29.2402460791314*(-0.184930902191118 - 0.446462692171690j +z)**2)
    elif n == 2:
        return abs(-29.2402460791314*(0.184930902191118 - 0.446462692171690j +z)**2)
    elif n == 3:
        return abs(-29.2402460791314*(0.446462692171690 - 0.184930902191118j +z)**2)
    elif n == 4:
        return abs(-29.2402460791314*(0.446462692171690 + 0.184930902191118j +z)**2)
    elif n == 5:
        return abs(-29.2402460791314*(0.184930902191118 + 0.446462692171690j +z)**2)
    elif n == 6:
        return abs(-29.2402460791314*(-0.184930902191118 + 0.446462692171690j + z)**2)
    elif n == 7:
        return abs(-29.2402460791314*(-0.446462692171690 + 0.184930902191118j +z)**2)
    elif n == 8:
        return abs(-5.8284271247462*(-1.00000000000000 - 0.414213562373095j + z)**2)
    elif n == 9:
        return abs(-5.8284271247462*(-0.414213562373095 - 1.00000000000000j + z)**2)
    elif n == 10:
        return abs(-5.8284271247462*(0.414213562373095 - 1.00000000000000j + z)**2)
    elif n == 11:
        return abs(-5.8284271247462*(1.00000000000000 - 0.414213562373095j + z)**2)
    elif n == 12:
        return abs(-5.8284271247462*(1.00000000000000 + 0.414213562373095j + z)**2)
    elif n == 13:
        return abs(-5.8284271247462*(0.414213562373095 + 1.00000000000000j + z)**2)
    elif n == 14:
        return abs(-5.8284271247462*(-0.414213562373095 + 1.00000000000000j + z)**2)
    elif n == 15:
        return abs(-5.8284271247462*(-1.00000000000000 + 0.414213562373095j + z)**2)

def samplePoint(word):
    points = [[0.468786, 0.194177], [0.194177, 0.468786], [-0.194177, 
  0.468786], [-0.468786, 
  0.194177], [-0.468786, -0.194177], [-0.194177, -0.468786], \
[0.194177, -0.468786], [0.468786, -0.194177], [0.8, 
  0.331371], [0.331371, 0.8], [-0.331371, 0.8], [-0.8, 
  0.331371], [-0.8, -0.331371], [-0.331371, -0.8], [0.331371, -0.8], \
[0.8, -0.331371]]
    p = points[word[-1]][0] + points[word[-1]][1]*1j
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    d = derivatives(word[0], samplePoint(word))
    if d>1:
        print(d, word)
    return d

def generateTree(words, dc):
    tt = time.time()
    generators = [np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 
  0., 0., 0., 0.], [2.70711, 5.41421, 0., 0., -0.292893, 0., 0., 0., 
  6.82843, 0.], [7.53553, 10.6569, 0., 0., -0.707107, 0., 0., 0., 
  16.4853, 0.], [11.6569, 13.6569, 0., 0., -1., 0., 0., 0., 23.3137, 
  0.], [12.6569, 12.6569, 0., 0., -1., 0., 0., 0., 23.3137, 
  0.], [9.94975, 8.24264, 0., 0., -0.707107, 0., 0., 0., 16.4853, 
  0.], [5.12132, 3., 0., 0., -0.292893, 0., 0., 0., 6.82843, 0.], [0.,
   0., 0., 0., 0., 0., 0., 0., 1., 0.], [2.17157, 2.34315, 0., 
  0., -0.171573, 0., 0., 0., 3., 0.]]),np.array([[0., 4.82843, 3.41421, 
  0., -0.414214, 0., 0., 0., 6.82843, 0.], [0., 
   1., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0., 
   0., 0., 0.], [0., 2.41421, 5.82843, 0., -0.414214, 0., 0., 0., 
   6.82843, 0.], [0., 6.82843, 11.6569, 0., -1., 0., 0., 0., 16.4853, 
   0.], [0., 10.6569, 15.0711, 0., -1.41421, 0., 0., 0., 23.3137, 
   0.], [0., 11.6569, 14.0711, 0., -1.41421, 0., 0., 0., 23.3137, 
   0.], [0., 9.24264, 9.24264, 0., -1., 0., 0., 0., 16.4853, 0.], [0.,
    0., 0., 0., 0., 0., 0., 0., 1., 0.], [0., 2., 2.58579, 
   0., -0.242641, 0., 0., 0., 3., 0.]]), np.array([[-1., 0., 11.6569, 6.82843, 
   0., 0., 0., 0., 16.4853, 0.], [-0.414214, 0., 5.82843, 2.41421, 0.,
    0., 0., 0., 6.82843, 0.], [0., 0., 1., 0., 0., 0., 0., 0., 0., 
   0.], [0., 0., 0., 1., 0., 0., 0., 0., 0., 0.], [-0.414214, 0., 
   3.41421, 4.82843, 0., 0., 0., 0., 6.82843, 0.], [-1., 0., 9.24264, 
   9.24264, 0., 0., 0., 0., 16.4853, 0.], [-1.41421, 0., 14.0711, 
   11.6569, 0., 0., 0., 0., 23.3137, 0.], [-1.41421, 0., 15.0711, 
   10.6569, 0., 0., 0., 0., 23.3137, 0.], [0., 0., 0., 0., 0., 0., 0.,
    0., 1., 0.], [-0.242641, 0., 2.58579, 2., 0., 0., 0., 0., 3., 
   0.]]), np.array([[0., 0., 0., 11.6569, 14.0711, 0., -1.41421, 0., 23.3137, 
   0.], [0., 0., 0., 9.24264, 9.24264, 0., -1., 0., 16.4853, 0.], [0.,
    0., 0., 4.82843, 3.41421, 0., -0.414214, 0., 6.82843, 0.], [0., 
   0., 0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 1., 0., 0., 
   0., 0., 0.], [0., 0., 0., 2.41421, 5.82843, 0., -0.414214, 0., 
   6.82843, 0.], [0., 0., 0., 6.82843, 11.6569, 0., -1., 0., 16.4853, 
   0.], [0., 0., 0., 10.6569, 15.0711, 0., -1.41421, 0., 23.3137, 
   0.], [0., 0., 0., 0., 0., 0., 0., 0., 1., 0.], [0., 0., 0., 2., 
   2.58579, 0., -0.242641, 0., 3., 0.]]), np.array([[0., 0., 0., -3.41421, 
   17.4853, 10.2426, 0., 0., 23.3137, 0.], [0., 0., 0., -3.41421, 
   18.4853, 9.24264, 0., 0., 23.3137, 0.], [0., 0., 0., -2.41421, 
   14.0711, 5.82843, 0., 0., 16.4853, 0.], [0., 0., 0., -1., 6.82843, 
   2., 0., 0., 6.82843, 0.], [0., 0., 0., 0., 1., 0., 0., 0., 0., 
   0.], [0., 0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., -1., 
   4.41421, 4.41421, 0., 0., 6.82843, 0.], [0., 0., 0., -2.41421, 
   11.6569, 8.24264, 0., 0., 16.4853, 0.], [0., 0., 0., 0., 0., 0., 
   0., 0., 1., 0.], [0., 0., 0., -0.585786, 3.17157, 1.75736, 0., 0., 
   3., 0.]]), np.array([[0., 0., 0., 0., -2.41421, 11.6569, 8.24264, 0., 
   16.4853, 0.], [0., 0., 0., 0., -3.41421, 17.4853, 10.2426, 0., 
   23.3137, 0.], [0., 0., 0., 0., -3.41421, 18.4853, 9.24264, 0., 
   23.3137, 0.], [0., 0., 0., 0., -2.41421, 14.0711, 5.82843, 0., 
   16.4853, 0.], [0., 0., 0., 0., -1., 6.82843, 2., 0., 6.82843, 
   0.], [0., 0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 0., 0.,
    0., 1., 0., 0., 0.], [0., 0., 0., 0., -1., 4.41421, 4.41421, 0., 
   6.82843, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 1., 0.], [0., 0., 
   0., 0., -0.585786, 3.17157, 1.75736, 0., 3., 0.]]), np.array([[0., 0., 
   0., -0.292893, 0., 0., 3., 5.12132, 6.82843, 0.], [0., 0., 
   0., -0.707107, 0., 0., 8.24264, 9.94975, 16.4853, 0.], [0., 0., 
   0., -1., 0., 0., 12.6569, 12.6569, 23.3137, 0.], [0., 0., 0., -1., 
   0., 0., 13.6569, 11.6569, 23.3137, 0.], [0., 0., 0., -0.707107, 0.,
    0., 10.6569, 7.53553, 16.4853, 0.], [0., 0., 0., -0.292893, 0., 
   0., 5.41421, 2.70711, 6.82843, 0.], [0., 0., 0., 0., 0., 0., 1., 
   0., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 1., 0., 0.], [0., 0., 0.,
    0., 0., 0., 0., 0., 1., 0.], [0., 0., 0., -0.171573, 0., 0., 
   2.34315, 2.17157, 3., 0.]]),np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [5.82843, 0., -0.414214, 
  0., 0., 0., 0., 2.41421, 6.82843, 0.], [11.6569, 0., -1., 0., 0., 
  0., 0., 6.82843, 16.4853, 0.], [15.0711, 0., -1.41421, 0., 0., 0., 
  0., 10.6569, 23.3137, 0.], [14.0711, 0., -1.41421, 0., 0., 0., 0., 
  11.6569, 23.3137, 0.], [9.24264, 0., -1., 0., 0., 0., 0., 9.24264, 
  16.4853, 0.], [3.41421, 0., -0.414214, 0., 0., 0., 0., 4.82843, 
  6.82843, 0.], [0., 0., 0., 0., 0., 0., 0., 1., 0., 0.], [0., 0., 0.,
   0., 0., 0., 0., 0., 1., 0.], [2.58579, 0., -0.242641, 0., 0., 0., 
  0., 2., 3., 0.]]),np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [5.41421, 0., 0., 0., 0., 
  0., 0., 3., -1.70711, 5.12132], [10.6569, 0., 0., 0., 0., 0., 0., 
  8.24264, -4.12132, 12.364], [13.6569, 0., 0., 0., 0., 0., 0., 
  12.6569, -5.82843, 17.4853], [12.6569, 0., 0., 0., 0., 0., 0., 
  13.6569, -5.82843, 17.4853], [8.24264, 0., 0., 0., 0., 0., 0., 
  10.6569, -4.12132, 12.364], [3., 0., 0., 0., 0., 0., 0., 
  5.41421, -1.70711, 5.12132], [0., 0., 0., 0., 0., 0., 0., 1., 0., 
  0.], [2.34315, 0., 0., 0., 0., 0., 0., 2.34315, -1., 2.], [0., 0., 
  0., 0., 0., 0., 0., 0., 0., 1.]]),np.array([[0., 5.82843, 2.41421, 0., 0., 0., 0., -0.414214, 0., 6.82843], [0., 
   1., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0., 
   0., 0., 0.], [0., 3.41421, 4.82843, 0., 0., 0., 0., -0.414214, 0., 
   6.82843], [0., 9.24264, 9.24264, 0., 0., 0., 0., -1., 0., 
   16.4853], [0., 14.0711, 11.6569, 0., 0., 0., 0., -1.41421, 0., 
   23.3137], [0., 15.0711, 10.6569, 0., 0., 0., 0., -1.41421, 0., 
   23.3137], [0., 11.6569, 6.82843, 0., 0., 0., 0., -1., 0., 
   16.4853], [0., 2.58579, 2., 0., 0., 0., 0., -0.242641, 0., 
   3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]), np.array([[-1., 0., 11.6569,
    6.82843, 0., 0., 0., 0., 0., 16.4853], [-0.414214, 0., 5.82843, 
   2.41421, 0., 0., 0., 0., 0., 6.82843], [0., 0., 1., 0., 0., 0., 0.,
    0., 0., 0.], [0., 0., 0., 1., 0., 0., 0., 0., 0., 0.], [-0.414214,
    0., 3.41421, 4.82843, 0., 0., 0., 0., 0., 6.82843], [-1., 0., 
   9.24264, 9.24264, 0., 0., 0., 0., 0., 16.4853], [-1.41421, 0., 
   14.0711, 11.6569, 0., 0., 0., 0., 0., 23.3137], [-1.41421, 0., 
   15.0711, 10.6569, 0., 0., 0., 0., 0., 23.3137], [-0.242641, 0., 
   2.58579, 2., 0., 0., 0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 
   0., 0., 1.]]), np.array([[0., 0., -3.41421, 18.4853, 9.24264, 0., 0., 0., 0.,
    23.3137], [0., 0., -2.41421, 14.0711, 5.82843, 0., 0., 0., 0., 
   16.4853], [0., 0., -1., 6.82843, 2., 0., 0., 0., 0., 6.82843], [0.,
    0., 0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 1., 0., 0., 
   0., 0., 0.], [0., 0., -1., 4.41421, 4.41421, 0., 0., 0., 0., 
   6.82843], [0., 0., -2.41421, 11.6569, 8.24264, 0., 0., 0., 0., 
   16.4853], [0., 0., -3.41421, 17.4853, 10.2426, 0., 0., 0., 0., 
   23.3137], [0., 0., -0.585786, 3.17157, 1.75736, 0., 0., 0., 0., 
   3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]), np.array([[-1., 0., 0., 0., 
   11.6569, 13.6569, 0., 0., 0., 23.3137], [-1., 0., 0., 0., 12.6569, 
   12.6569, 0., 0., 0., 23.3137], [-0.707107, 0., 0., 0., 9.94975, 
   8.24264, 0., 0., 0., 16.4853], [-0.292893, 0., 0., 0., 5.12132, 3.,
    0., 0., 0., 6.82843], [0., 0., 0., 0., 1., 0., 0., 0., 0., 
   0.], [0., 0., 0., 0., 0., 1., 0., 0., 0., 0.], [-0.292893, 0., 0., 
   0., 2.70711, 5.41421, 0., 0., 0., 6.82843], [-0.707107, 0., 0., 0.,
    7.53553, 10.6569, 0., 0., 0., 16.4853], [-0.171573, 0., 0., 0., 
   2.17157, 2.34315, 0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0.,
    0., 1.]]), np.array([[0., 0., 0., 0., 0., 5.82843, 14.0711, -2.41421, 0., 
   16.4853], [0., 0., 0., 0., 0., 9.24264, 18.4853, -3.41421, 0., 
   23.3137], [0., 0., 0., 0., 0., 10.2426, 17.4853, -3.41421, 0., 
   23.3137], [0., 0., 0., 0., 0., 8.24264, 11.6569, -2.41421, 0., 
   16.4853], [0., 0., 0., 0., 0., 4.41421, 4.41421, -1., 0., 
   6.82843], [0., 0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 
   0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 0., 0., 2., 6.82843, -1.,
    0., 6.82843], [0., 0., 0., 0., 0., 1.75736, 3.17157, -0.585786, 
   0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 
   1.]]), np.array([[0., -0.414214, 0., 0., 0., 0., 2.41421, 5.82843, 0., 
   6.82843], [0., -1., 0., 0., 0., 0., 6.82843, 11.6569, 0., 
   16.4853], [0., -1.41421, 0., 0., 0., 0., 10.6569, 15.0711, 0., 
   23.3137], [0., -1.41421, 0., 0., 0., 0., 11.6569, 14.0711, 0., 
   23.3137], [0., -1., 0., 0., 0., 0., 9.24264, 9.24264, 0., 
   16.4853], [0., -0.414214, 0., 0., 0., 0., 4.82843, 3.41421, 0., 
   6.82843], [0., 0., 0., 0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 
   0., 0., 0., 0., 1., 0., 0.], [0., -0.242641, 0., 0., 0., 0., 2., 
   2.58579, 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
   np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 
  0., 0., 0., 0.], [2., 6.82843, -1., 0., 0., 0., 0., 0., 0., 
  6.82843], [5.82843, 14.0711, -2.41421, 0., 0., 0., 0., 0., 0., 
  16.4853], [9.24264, 18.4853, -3.41421, 0., 0., 0., 0., 0., 0., 
  23.3137], [10.2426, 17.4853, -3.41421, 0., 0., 0., 0., 0., 0., 
  23.3137], [8.24264, 11.6569, -2.41421, 0., 0., 0., 0., 0., 0., 
  16.4853], [4.41421, 4.41421, -1., 0., 0., 0., 0., 0., 0., 
  6.82843], [1.75736, 3.17157, -0.585786, 0., 0., 0., 0., 0., 0., 
  3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]])]

    root = Node([3.61312592975275, 3.61312592975275, 3.61312592975275, \
3.61312592975275, 3.61312592975275, 3.61312592975275, \
3.61312592975275, 3.61312592975275, 2.23982880884355, \
-1.00000000000000], [], words, False)
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
    matrix = matrix.power(a)
    vec = np.ones(l)
    previous_entry = vec[0]
    previous_val = 0
    current = matrix * vec
    current_val = current[0] / previous_entry
    count = 0
    while count < 10000000000 and abs(current_val - previous_val) > 1e-15:
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
    matrix = constructMatrix(words,m)
    #print(matrix)

    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*15)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()