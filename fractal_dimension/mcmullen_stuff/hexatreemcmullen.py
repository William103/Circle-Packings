from hexanode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix
from scipy import stats

def functions(n, z):
    if n == 0:
        return 0.333333333333333 + 0.192450089729875j + np.conj(0.0370370370370370/(-0.333333333333333 - 0.192450089729875j + z))
    elif n == 1:
        return 0.384900179459751j + np.conj(0.0370370370370370/(-0.384900179459751j + z))
    elif n == 2:
        return -0.333333333333333 + 0.192450089729875j + np.conj(0.0370370370370370/(0.333333333333333 - 0.192450089729875j + z))
    elif n == 3:
        return -0.333333333333333 - 0.192450089729875j + np.conj(0.0370370370370370/(0.333333333333333 + 0.192450089729875j + z))
    elif n == 4:
        return -0.384900179459751j + np.conj(0.0370370370370370/(0.384900179459751j + z))
    elif n == 5:
        return 0.333333333333333 - 0.192450089729875j + np.conj(0.0370370370370370/(-0.333333333333333 + 0.192450089729875j + z))
    elif n == 6:
        return 1.00000000000000 + 0.577350269189626j + np.conj(0.333333333333333/(-1.00000000000000 - 0.577350269189626j + z))
    elif n == 7:
        return 1.15470053837925j + np.conj(0.333333333333333/(-1.15470053837925j + z))
    elif n == 8:
        return -1.00000000000000 + 0.577350269189626j + np.conj(0.333333333333333/(1.00000000000000 - 0.577350269189626j + z))
    elif n == 9:
        return -1.00000000000000 - 0.577350269189626j + np.conj(0.333333333333333/(1.00000000000000 + 0.577350269189626j + z))
    elif n == 10:
        return -1.15470053837925j + np.conj(0.333333333333333/(1.15470053837925j + z))
    elif n == 11:
        return 1.00000000000000 - 0.577350269189626j + np.conj(0.333333333333333/(-1.00000000000000 + 0.577350269189626j + z))

def derivatives(n, z):
    if n == 0:
        return abs(-27.0000000000000*(-0.333333333333333 - 0.192450089729875j + z)**2)
    elif n == 1:
        return abs(-27.0000000000000*(-0.384900179459751j + z)**2)
    elif n == 2:
        return abs(-27.0000000000000*(0.333333333333333 - 0.192450089729875j + z)**2)
    elif n == 3:
        return abs(-27.0000000000000*(0.333333333333333 + 0.192450089729875j + z)**2)
    elif n == 4:
        return abs(-27.0000000000000*(0.384900179459751j + z)**2)
    elif n == 5:
        return abs(-27.0000000000000*(-0.333333333333333 + 0.192450089729875j + z)**2)
    elif n == 6:
        return abs(-3.00000000000000*(-1.00000000000000 - 0.577350269189626j + z)**2)
    elif n == 7:
        return abs(-3.00000000000000*(-1.15470053837925j + z)**2)
    elif n == 8:
        return abs(-3.00000000000000*(1.00000000000000 - 0.577350269189626j + z)**2)
    elif n == 9:
        return abs(-3.00000000000000*(1.00000000000000 + 0.577350269189626j + z)**2)
    elif n == 10:
        return abs(-3.00000000000000*(1.15470053837925j + z)**2)
    elif n == 11:
        return abs(-3.00000000000000*(-1.00000000000000 + 0.577350269189626j + z)**2)

def samplePoint(word):
    points = [0.333333333333333 + 0.192450089729875j, 
    0.384900179459751j, 
    -0.333333333333333 + 0.192450089729875j, 
    -0.333333333333333 - 0.192450089729875j, 
    -0.384900179459751j, 
    0.333333333333333 - 0.192450089729875j, 
    0.750000000000000 + 0.433012701892219j, 
    0.866025403784439j, 
    -0.750000000000000 + 0.433012701892219j, 
    -0.750000000000000 - 0.433012701892219j, 
    -0.866025403784439j, 
    0.750000000000000 - 0.433012701892219j]
    p = points[word[-1]]
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    generators = [np.array([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [2, 6, -1, 0, 0, 
  0, 6, 0], [5, 10, -2, 0, 0, 0, 12, 0], [6, 9, -2, 0, 0, 0, 12, 
  0], [4, 4, -1, 0, 0, 0, 6, 0], [0, 0, 0, 0, 0, 0, 1, 0], [2, 10/
  3, -(2/3), 0, 0, 0, 3, 0]]),
    np.array([[0, 5, 3, 0, 0, 0, 9/2, -(3/2)], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 
  0, 0, 0, 0, 0], [0, 3, 5, 0, 0, 0, 9/2, -(3/2)], [0, 7, 8, 0, 0, 0, 
  9, -3], [0, 8, 7, 0, 0, 0, 9, -3], [0, 0, 0, 0, 0, 0, 1, 0], [0, 8/
  3, 8/3, 0, 0, 0, 2, -1]]),
    np.array([[0, -2, 10, 5, 0, 0, 12, 0], [0, -1, 6, 2, 0, 0, 6, 0], [0, 0, 1, 0, 
  0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, -1, 4, 4, 0, 0, 6, 
  0], [0, -2, 9, 6, 0, 0, 12, 0], [0, 0, 0, 0, 0, 0, 1, 
  0], [0, -(2/3), 10/3, 2, 0, 0, 3, 0]]),
    np.array([[-1, 0, 0, 6, 8, 0, 12, 0], [-1, 0, 0, 7, 7, 0, 12, 0], [-(1/2), 0, 
  0, 9/2, 3, 0, 6, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0,
   0], [-(1/2), 0, 0, 5/2, 5, 0, 6, 0], [0, 0, 0, 0, 0, 0, 1, 
  0], [-(1/3), 0, 0, 7/3, 8/3, 0, 3, 0]]),
    np.array([[-1, 0, 0, 0, 2, 6, 6, 0], [-2, 0, 0, 0, 5, 10, 12, 0], [-2, 0, 0, 0,
   6, 9, 12, 0], [-1, 0, 0, 0, 4, 4, 6, 0], [0, 0, 0, 0, 1, 0, 0, 
  0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [-(2/3), 0, 
  0, 0, 2, 10/3, 3, 0]]),
    np.array([[1, 0, 0, 0, 0, 0, 0, 0], [4, 0, 0, 0, -1, 4, 6, 0], [6, 0, 0, 0, -2,
   9, 12, 0], [5, 0, 0, 0, -2, 10, 12, 0], [2, 0, 0, 0, -1, 6, 6, 
  0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [2, 0, 0, 
  0, -(2/3), 10/3, 3, 0]]),
    np.array([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0], [4, 4, 0, 0, 
  0, -1, 0, 6], [9, 6, 0, 0, 0, -2, 0, 12], [10, 5, 0, 0, 0, -2, 0, 
  12], [6, 2, 0, 0, 0, -1, 0, 6], [10/3, 2, 0, 0, 0, -(2/3), 0, 
  3], [0, 0, 0, 0, 0, 0, 0, 1]]),
    np.array([[0, 4, 4, -1, 0, 0, 0, 6], [0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 
  0, 0, 0], [0, 2, 6, -1, 0, 0, 0, 6], [0, 5, 10, -2, 0, 0, 0, 
  12], [0, 6, 9, -2, 0, 0, 0, 12], [0, 2, 10/3, -(2/3), 0, 0, 0, 
  3], [0, 0, 0, 0, 0, 0, 0, 1]]),
    np.array([[0, 0, 7, 7, 0, -1, 0, 12], [0, 0, 9/2, 3, 0, -(1/2), 0, 6], [0, 0, 
  1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 5/2, 5, 
  0, -(1/2), 0, 6], [0, 0, 6, 8, 0, -1, 0, 12], [0, 0, 7/3, 8/3, 
  0, -(1/3), 0, 3], [0, 0, 0, 0, 0, 0, 0, 1]]),
    np.array([[0, -1, 0, 7, 7, 0, 0, 12], [0, -1, 0, 8, 6, 0, 0, 12], [0, -(1/2), 
  0, 5, 5/2, 0, 0, 6], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0,
   0], [0, -(1/2), 0, 3, 9/2, 0, 0, 6], [0, -(1/3), 0, 8/3, 7/3, 0, 0,
   3], [0, 0, 0, 0, 0, 0, 0, 1]]),
    np.array([[0, 0, -(1/2), 0, 3, 9/2, 0, 6], [0, 0, -1, 0, 7, 7, 0, 12], [0, 
   0, -1, 0, 8, 6, 0, 12], [0, 0, -(1/2), 0, 5, 5/2, 0, 6], [0, 0, 0, 
   0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, -(1/3), 0, 8/3, 7/
   3, 0, 3], [0, 0, 0, 0, 0, 0, 0, 1]]),
    np.array([[1, 0, 0, 0, 0, 0, 0, 0], [5, 0, -(1/2), 0, 0, 5/2, 0, 6], [8, 0, -1,
   0, 0, 6, 0, 12], [7, 0, -1, 0, 0, 7, 0, 12], [3, 0, -(1/2), 0, 0, 
  9/2, 0, 6], [0, 0, 0, 0, 0, 1, 0, 0], [8/3, 0, -(1/3), 0, 0, 7/3, 0,
   3], [0, 0, 0, 0, 0, 0, 0, 1]])]

    root = Node([3, 3, 3, 3, 3, 3, 3, -1], [], words, False)
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

    m = 1000
    print("maximum:",m)
    
    matrix = constructMatrix(words,m)
    #print(matrix)
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*11)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()