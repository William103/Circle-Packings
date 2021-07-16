from pentanode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix
from scipy import stats

def functions(n, z):
    if n == 0:
        return 0.259616 + 0.188622j + np.conj(0.0355783/(-0.259616 - 0.188622j + z))
    elif n == 1:
        return -0.0991646 + 0.305197j + np.conj(0.0355783/(0.0991646 - 0.305197j + z))
    elif n == 2:
        return -0.320903 + np.conj(0.0355783/(0.320903 + z))
    elif n == 3:
        return -0.0991646 - 0.305197j + np.conj(0.0355783/(0.0991646 + 0.305197j + z))
    elif n == 4:
        return 0.259616 - 0.188622j + np.conj(0.0355783/(-0.259616 + 0.188622j + z))
    elif n == 5:
        return 1. + 0.726543j + np.conj(0.527864/(-1. - 0.726543j + z))
    elif n == 6:
        return -0.381966 + 1.17557j + np.conj(0.527864/(0.381966 - 1.17557j + z))
    elif n == 7:
        return -1.23607 + np.conj(0.527864/(1.23607 + z))
    elif n == 8:
        return -0.381966 - 1.17557j + np.conj(0.527864/(0.381966 + 1.17557j + z))
    elif n == 9:
        return 1. - 0.726543j + np.conj(0.527864/(-1. + 0.726543j + z))

def derivatives(n, z):
    if n == 0:
        return abs(-28.107*(-0.259616 - 0.188622j + z)**2)
    elif n == 1:
        return abs(-28.107*(0.0991646 - 0.305197j + z)**2)
    elif n == 2:
        return abs(-28.107*(0.320903 + z)**2)
    elif n == 3:
        return abs(-28.107*(0.0991646 + 0.305197j + z)**2)
    elif n == 4:
        return abs(-28.107*(-0.259616 + 0.188622j + z)**2)
    elif n == 5:
        return abs(-1.89443*(-1. - 0.726543j + z)**2)
    elif n == 6:
        return abs(-1.89443*(0.381966 - 1.17557j + z)**2)
    elif n == 7:
        return abs(-1.89443*(1.23607 + z)**2)
    elif n == 8:
        return abs(-1.89443*(0.381966 + 1.17557j + z)**2)
    elif n == 9:
        return abs(-1.89443*(-1. + 0.726543j + z)**2)

def samplePoint(word):
    if word[-1] == 0:
        p = 0.259616 + 0.188622j
    elif word[-1] == 1:
        p = -0.0991646 + 0.305197j
    elif word[-1] == 2:
        p = -0.320903
    elif word[-1] == 3:
        p = -0.0991646 - 0.305197j
    elif word[-1] == 4:
        p = 0.259616 - 0.188622j
    elif word[-1] == 5:
        p = 0.6 + 0.435926j
    elif word[-1] == 6:
        p = -0.22918 + 0.705342j
    elif word[-1] == 7:
        p = -0.741641
    elif word[-1] == 8:
        p = -0.22918 - 0.705342j
    elif word[-1] == 9:
        p = 0.6 - 0.435926j
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    generators = [np.array([[1., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0.], [2.61803, 4.23607, 0., -0.618034, 0., 5.23607, 0.], [5.23607, 5.23607, 0., -1., 0., 8.47214, 0.], [4.23607, 2.61803, 0., -0.618034, 0., 5.23607, 0.], [0., 0., 0., 0., 0., 1., 0.], [2.76393, 2.76393, 0., -0.472136, 0., 3., 0.]]),
    np.array([[0., 4.23607, 2.61803, 0., -0.618034, 5.23607, 0.], [0., 1., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0.], [0., 2.61803, 4.23607, 0., -0.618034, 5.23607, 0.], [0., 5.23607, 5.23607, 0., -1., 8.47214, 0.], [0., 0., 0., 0., 0., 1., 0.], [0., 2.76393, 2.76393, 0., -0.472136, 3., 0.]]),
    np.array([[0., 0., 4.23607, 6.8541, -1.61803, 8.47214, 0.], [0., 0., 3.61803, 3.61803, -1., 5.23607, 0.], [0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 1., 0., 0., 0.], [0., 0., 2., 5.23607, -1., 5.23607, 0.], [0., 0., 0., 0., 0., 1., 0.], [0., 0., 2.2918, 3.52786, -0.763932, 3., 0.]]),
    np.array([[0., -0.618034, 0., 2.61803, 4.23607, 5.23607, 0.], [0., -1., 0., 5.23607, 5.23607, 8.47214, 0.], [0., -0.618034, 0., 4.23607, 2.61803, 5.23607, 0.], [0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 0., 1., 0.], [0., -0.472136, 0., 2.76393, 2.76393, 3., 0.]]),
    np.array([[1., 0., 0., 0., 0., 0., 0.], [4.23607, 0., -0.618034, 0., 2.61803, 5.23607, 0.], [5.23607, 0., -1., 0., 5.23607, 8.47214, 0.], [2.61803, 0., -0.618034, 0., 4.23607, 5.23607, 0.], [0., 0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 0., 1., 0.], [2.76393, 0., -0.472136, 0., 2.76393, 3., 0.]]),
    np.array([[1., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0.], [2.61803,4.23607, 0., -0.618034, 0., 0., 5.23607], [5.23607, 5.23607, 0., -1., 0., 0., 8.47214], [4.23607, 2.61803, 0., -0.618034, 0., 0.,5.23607], [2.76393, 2.76393, 0., -0.472136, 0., 0., 3.], [0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[-1., 5.23607, 2., 0., 0., 0., 5.23607], [0., 1., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0.], [-1., 3.61803, 3.61803, 0., 0., 0., 5.23607], [-1.61803, 6.8541, 4.23607, 0., 0., 0., 8.47214], [-0.763932, 3.52786, 2.2918, 0., 0., 0., 3.], [0., 0., 0.,0., 0., 0., 1.]]),
    np.array([[-1., 0., 5.23607, 5.23607, 0., 0., 8.47214], [-0.618034, 0., 4.23607, 2.61803, 0., 0., 5.23607], [0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 1., 0., 0., 0.], [-0.618034, 0., 2.61803, 4.23607,0., 0., 5.23607], [-0.472136, 0., 2.76393, 2.76393, 0., 0., 3.], [0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[-1., 0., 0., 2., 5.23607, 0., 5.23607], [-1.61803, 0., 0., 4.23607, 6.8541, 0., 8.47214], [-1., 0., 0., 3.61803, 3.61803, 0., 5.23607], [0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 0., 1., 0., 0.], [-0.763932, 0., 0., 2.2918, 3.52786, 0., 3.], [0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[1., 0., 0., 0., 0., 0., 0.], [4.61803, 0., 0., 0., 3., -1.30902, 3.92705], [5.8541, 0., 0., 0., 5.8541, -2.11803, 6.3541], [3., 0., 0., 0., 4.61803, -1.30902, 3.92705], [0., 0., 0., 0., 1., 0., 0.], [3.05573, 0., 0., 0., 3.05573, -1., 2.], [0., 0., 0., 0., 0., 0., 1.]])]
    root = Node([2.7013, 2.7013, 2.7013, 2.7013, 2.7013, 3.85184, -1.], [], words, False)
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
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*9)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()