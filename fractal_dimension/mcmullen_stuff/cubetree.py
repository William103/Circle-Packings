#https://www.desmos.com/calculator/vaskxmhoxh

from cubenode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix

def functions(n, z):
    if n == 0:
        return -1.732050808 + 1j + np.conj(3./(1.732050808 - 1.000000000j + z))
    elif n == 1:
        return -2.000000000j + np.conj(3./(2.000000000j + z))
    elif n == 2:
        return 1.732050808 + 1j + np.conj(3./(-1.732050808 - 1.000000000j + z))
    elif n == 3:
        return 0.202041j + np.conj(0.0306154/(-0.202041j + z))
    elif n == 4:
        return 0.174973 - 0.101021j + np.conj(0.0306154/(-0.174973 + 0.101021j + z))
    elif n == 5:
        return -0.174973 - 0.101021j + np.conj(0.0306154/(0.174973 + 0.101021j + z))

def derivatives(n, z):
    if n == 0:
        return abs(-0.333333*(1.732050808 - 1.000000000j + z)**2)
    elif n == 1:
        return abs(-0.333333*(2.000000000j + z)**2)
    elif n == 2:
        return abs(-0.333333*(-1.732050808 - 1.000000000j + z)**2)
    elif n == 3:
        return abs(-32.6633*(-0.202041j + z)**2)
    elif n == 4:
        return abs(-32.6633*(-0.174973 + 0.101021j + z)**2)
    elif n == 5:
        return abs(-32.6633*(0.174973 + 0.101021j + z)**2)

def samplePoint(word):
    points = [-0.606218 + 0.35j, 0. - 0.7j, 0.606218 + 0.35j, 0. + 
 0.181837j, 0.157475 - 0.0909185j, -0.157475 - 0.0909185j]
    p = points[word[-1]]
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    
    generators = [np.array([[0., 1., 1., 0., -1., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0., 
  0.], [0., 0., 1., 0., 0., 0., 0., 0.], [0., 3., 3., 0., 0., 0., 
  0., -1.], [0., 0., 0., 0., 1., 0., 0., 0.], [0., 3., 2., 0., 1., 0.,
   0., -1.], [0., 2., 3., 0., 1., 0., 0., -1.], [0., 2., 2., 0., 2., 
  0., 0., -1.]]), 
  np.array([[1., 0., 0., 0., 0., 0., 0., 0.], [1., 0., 3., 2., -1., 0., 0., 
  0.], [0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 0., 1., 0., 0., 0., 
  0.], [0., 0., 4., 2., -1., 0., 0., 0.], [0., 0., 3., 3., -1., 0., 
  0., 0.], [-1., 0., 1., 1., 0., 0., 0., 0.], [-1., 0., 4., 3., -1., 
  0., 0., 0.]]),
  np.array([[1., 0., 0., 0., 0., 0., 0., 0.], [1., 0., 0., -1., 0., 1., 0., 
  0.], [4., 0., 0., -1., -1., 3., 0., 0.], [0., 0., 0., 1., 0., 0., 
  0., 0.], [4., 0., 0., -2., -1., 4., 0., 0.], [0., 0., 0., 0., 0., 
  1., 0., 0.], [3., 0., 0., 0., -1., 3., 0., 0.], [3., 0., 
  0., -1., -1., 4., 0., 0.]]),
  np.array([[0., 0., 0., 0., 3., 3., -1., 0.], [0., 1., 0., 0., 0., 0., 0., 
  0.], [0., -1., 0., 0., 4., 3., -1., 0.], [0., -1., 0., 0., 3., 
  4., -1., 0.], [0., 0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 0., 0.,
   1., 0., 0.], [0., -2., 0., 0., 4., 4., -1., 0.], [0., -1., 0., 0., 
  1., 1., 0., 0.]]),
  np.array([[0., 0., 0., 0., -1., 3., 3., 0.], [0., 0., 0., -1., -1., 4., 3., 
  0.], [0., 0., 0., -1., -1., 3., 4., 0.], [0., 0., 0., 1., 0., 0., 
  0., 0.], [0., 0., 0., -2., -1., 4., 4., 0.], [0., 0., 0., 0., 0., 
  1., 0., 0.], [0., 0., 0., 0., 0., 0., 1., 0.], [0., 0., 0., -1., 0.,
  1., 1., 0.]]),
  np.array([[0., -1., 0., 0., 4., 0., 3., -1.], [0., -1., 0., 0., 4., 0., 2., 
  0.], [0., 0., 0., 0., 1., 0., 1., -1.], [0., -1., 0., 0., 3., 0., 
  3., 0.], [0., 0., 0., 0., 1., 0., 0., 0.], [0., -1., 0., 0., 3., 0.,
  2., 1.], [0., 0., 0., 0., 0., 0., 1., 0.], [0., 0., 0., 0., 0., 0.,
  0., 1.]])]

    root = Node([-1,2.632993161855454,2.632993161855454,2.632993161855454,6.265986323710909,6.265986323710909,6.265986323710909,9.898979485566363], [], words, False)
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
        count += 1
        print(count,k1,y1)

def main():
    start = time.time()

    m = 10
    print("maximum:",m)
    words = {}
    matrix = constructMatrix(words,m)
    #print(matrix)
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*5)
    secantMethod(matrix,matrix.shape[0],1,1.47,1.49,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()