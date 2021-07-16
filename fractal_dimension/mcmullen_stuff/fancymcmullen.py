from node import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix
from scipy import stats

def functions(n, z):
    if n == 0:
        return np.conj((2-math.sqrt(3))**2 / z)
    elif n ==1:
        return 2j + np.conj(3 / (z - 2j))
    elif n == 2:
        return -1 * math.sqrt(3) - 1j + np.conj(3 / (z + math.sqrt(3) + 1j))
    elif n == 3:
        return math.sqrt(3) - 1j + np.conj(3 / (z - math.sqrt(3) + 1j))

def derivatives(n, z):
    if n == 0:
        return abs(z**2 / (2-math.sqrt(3))**2)
    elif n == 1:
        return abs((z-2j)**2 / 3)
    elif n == 2:
        return abs((z+math.sqrt(3)+1j)**2 / 3)
    elif n == 3:
        return abs((z-math.sqrt(3)+1j)**2 / 3)

def samplePoint(word):
    if word[-1] == 0:
        p = 0 + 0j
    elif word[-1] == 1:
        p = 0 + 0.5j
    elif word[-1] == 2:
        p = -1 * math.sqrt(3) / 4 - 0.25j
    elif word[-1] == 3:
        p = math.sqrt(3) / 4 - 0.25j
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    generators = [np.array([[-1,2,2,2],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),
        np.array([[1,0,0,0],[2,-1,2,2],[0,0,1,0],[0,0,0,1]]),
        np.array([[1,0,0,0],[0,1,0,0],[2,2,-1,2],[0,0,0,1]]),
        np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[2,2,2,-1]])]
    root = Node([-1., 2. + math.sqrt(3), 2. + math.sqrt(3), 2. + math.sqrt(3)], [], words, False)
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
            for wor in thing.leaves():
                row.append(i)
                col.append(words[str(wor.word)])
                data.append(sampleValue(leaf.word))
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
    while count < 1000000000 and abs(current_val - previous_val) > 1e-15:
        previous_val = current_val
        previous_entry = current[0]
        current = matrix * current
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

    matrix = constructMatrix(words,400)
    #print(matrix.toarray())
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*3)
    secantMethod(matrix,matrix.shape[0],1,1.30,1.31,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()