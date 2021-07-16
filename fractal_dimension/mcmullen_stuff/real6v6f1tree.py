#https://www.desmos.com/calculator/1aoycwgvuz

from real6v6f1node import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix

def functions(n, z):
    if n == 0:
        return -1*np.conj(z)
    else:
        duals = [[],[0.60056621200156+1j,0.60056621200156],
        [1.11757463364131127+0.028662703870927425j,0.49979457015027892407],
        [0.63600982475703-1j,0.63600982475703],
        [0.312309249205871+0.133830541363598j,0.312309249205871],
        [0.174700306168316-0.333333333333333j,0.174700306168316]]
        return duals[n][0] + np.conj((duals[n][1])**2 / (z-duals[n][0]))

def derivatives(n, z):
    if n == 0:
        return abs(-1)
    else:
        duals = [[],[0.60056621200156+1j,0.60056621200156],
        [1.11757463364131127+0.028662703870927425j,0.49979457015027892407],
        [0.63600982475703-1j,0.63600982475703],
        [0.312309249205871+0.133830541363598j,0.312309249205871],
        [0.174700306168316-0.333333333333333j,0.174700306168316]]
        return abs(-1*((z-duals[n][0])**2) / ((duals[n][1])**2))

def samplePoint(word):
    points = [-0.5, 0.5+0.7j, 0.9+0.05j, 0.45-0.75j, 0.3+0.1j, 0.15-0.33j]
    p = points[word[-1]]
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    generators = [np.array([[0., -0.618034, 0., 2.23607, -0.527864, 2.8541], [0., -1., 0., 
  5.23607, -2.47214, 7.23607], [0., 0., 0., 1.61803, -1., 
  1.61803], [0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 1., 0.], [0., 
  0., 0., 0., 0., 1.]]),
  np.array([[0., 2.61803, 2.61803, 0., -0.381966, 3.8541], [0., 1., 0., 0., 0., 
  0.], [0., 0., 1., 0., 0., 0.], [0., 3.23607, 5.8541, 0., -0.618034, 
  6.23607], [0., 5.23607, 8.47214, 0., -1., 11.7082], [0., 0., 0., 0.,
   0., 1.]]),
   np.array([[1., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0.], [3., 2.61803, 
  0., -0.618034, 0., 3.23607], [6.47214, 3.23607, 0., -1., 0., 
  5.23607], [7.47214, 2.61803, 0., -1., 0., 6.8541], [0., 0., 0., 0., 
  0., 1.]]),
  np.array([[1., 0., 0., 0., 0., 0.], [6.8541, 0., 0., -1.61803, 3.61803, 
  5.8541], [5.23607, 0., 0., -1.61803, 4.23607, 4.8541], [3.23607, 0.,
   0., -1., 3.23607, 2.], [0., 0., 0., 0., 1., 0.], [0., 0., 0., 0., 
  0., 1.]]),
  np.array([[1., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0.], [-1., 0.618034, 
  0., 0.618034, 0., 0.], [0., 0., 0., 1., 0., 0.], [1., 1.38197, 0., 
  2.23607, 0., -1.61803], [0., 1.23607, 0., 0.763932, 0., -1.]]),
  np.array([[1., 0., 0., 0., 0., 0.], [10.0902, 0., 0., 6.8541, 
  6.8541, -2.61803], [5.23607, 0., 0., 4.8541, 
  4.23607, -1.61803], [0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 1., 
  0.], [3.23607, 0., 0., 2., 3.23607, -1.]])]
    root = Node([2.73606797749979, 3.88196601125011, 2.30901699437495, \
4.28115294937453, 3.00000000000000, -1.00000000000000], [], words, False)
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

    m = 500
    print("maximum:",m)
    matrix = constructMatrix(words,m)
    #print(matrix)
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*5)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()