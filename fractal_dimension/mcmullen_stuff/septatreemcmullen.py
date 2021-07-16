from septanode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix

def functions(n, z):
    if n == 0:
        return 0.394813223302777 + 0.190132027512207j + np.conj(0.0361501878859025/(-0.394813223302777 - 0.190132027512207j + z))
    elif n == 1:
        return 0.0975108134337358 + 0.427222787833377j + np.conj(0.0361501878859025/(-0.0975108134337358 - 0.427222787833377j + z))
    elif n == 2:
        return -0.273219227809010 + 0.342606075159329j + np.conj(0.0361501878859025/(0.273219227809010 - 0.342606075159329j + z))
    elif n == 3:
        return -0.438209617855007 + np.conj(0.0361501878859025/(0.438209617855007 + z))
    elif n == 4:
        return -0.273219227809010 - 0.342606075159329j + np.conj(0.0361501878859025/(0.273219227809010 + 0.342606075159329j + z))
    elif n == 5:
        return 0.0975108134337358 - 0.427222787833377j + np.conj(0.0361501878859025/(-0.0975108134337358 + 0.427222787833377j + z))
    elif n == 6:
        return 0.394813223302777 - 0.190132027512207j + np.conj(0.0361501878859025/(-0.394813223302777 + 0.190132027512207j + z))
    elif n == 7:
        return 1.00000000000000 + 0.481574618807529j + np.conj(0.231914113479617/(-1.00000000000000 - 0.481574618807529j + z))
    elif n == 8:
        return 0.246979603717467 + 1.08208834612853j + np.conj(0.231914113479617/(-0.246979603717467 - 1.08208834612853j + z))
    elif n == 9:
        return -0.692021471630096 + 0.867767478235116j + np.conj(0.231914113479617/(0.692021471630096 - 0.867767478235116j + z))
    elif n == 10:
        return -1.10991626417474 + np.conj(0.231914113479617/(1.10991626417474 + z))
    elif n == 11:
        return -0.692021471630096 - 0.867767478235116j + np.conj(0.231914113479617/(0.692021471630096 + 0.867767478235116j + z))
    elif n == 12:
        return 0.246979603717467 - 1.08208834612853j + np.conj(0.231914113479617/(-0.246979603717467 + 1.08208834612853j + z))
    elif n == 13:
        return 1.00000000000000 - 0.481574618807529j + np.conj(0.231914113479617/(-1.00000000000000 + 0.481574618807529j + z))

def derivatives(n, z):
    if n == 0:
        return abs(-27.6623735167354*(-0.394813223302777 - 0.190132027512207j + z)**2)
    elif n == 1:
        return abs(-27.6623735167354*(-0.0975108134337358 - 0.427222787833377j + z)**2)
    elif n == 2:
        return abs(-27.6623735167354*(0.273219227809010 - 0.342606075159329j + z)**2)
    elif n == 3:
        return abs(-27.6623735167354*(0.438209617855007 + z)**2)
    elif n == 4:
        return abs(-27.6623735167354*(0.273219227809010 + 0.342606075159329j + z)**2)
    elif n == 5:
        return abs(-27.6623735167354*(-0.0975108134337358 + 0.427222787833377j + z)**2)
    elif n == 6:
        return abs(-27.6623735167354*(-0.394813223302777 + 0.190132027512207j + z)**2)
    elif n == 7:
        return abs(-4.31194111042273*(-1.00000000000000 - 0.481574618807529j + z)**2)
    elif n == 8:
        return abs(-4.31194111042273*(-0.246979603717467 - 1.08208834612853j + z)**2)
    elif n == 9:
        return abs(-4.31194111042273*(0.692021471630096 - 0.867767478235116j + z)**2)
    elif n == 10:
        return abs(-4.31194111042273*(1.10991626417474 + z)**2)
    elif n == 11:
        return abs(-4.31194111042273*(0.692021471630096 + 0.867767478235116j + z)**2)
    elif n == 12:
        return abs(-4.31194111042273*(-0.246979603717467 + 1.08208834612853j + z)**2)
    elif n == 13:
        return abs(-4.31194111042273*(-1.00000000000000 + 0.481574618807529j + z)**2)


def samplePoint(word):
    points = [0.394813223302777 + 0.190132027512207j, 0.0975108134337358 + 
 0.427222787833377j, -0.273219227809010 + 
 0.342606075159329j, -0.438209617855007, -0.273219227809010 - 
 0.342606075159329j, 0.0975108134337358 - 
 0.427222787833377j, 0.394813223302777 - 0.190132027512207j, 0.8 + 
 0.38526j, 0.197584 + 0.865671j, -0.553617 + 
 0.694214j, -0.887933, -0.553617 - 0.694214j, 0.197584 - 
 0.865671j, 0.8 - 0.38526j]
    p = points[word[-1]]
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    
    generators = [np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0., 
  0., 0.], [3.24698, 4.69202, 0., 0., 0., -0.445042, 0., 6.49396, 
  0.], [8.2959, 8.2959, 0., 0., 0., -1., 0., 14.5918, 0.], [11.3448, 
  9.09783, 0., 0., 0., -1.24698, 0., 18.1957, 0.], [10.0978, 6.49396, 
  0., 0., 0., -1., 0., 14.5918, 0.], [5.49396, 2.44504, 0., 0., 
  0., -0.445042, 0., 6.49396, 0.], [0., 0., 0., 0., 0., 0., 0., 1., 
  0.], [2.61596, 2.122, 0., 0., 0., -0.274127, 0., 3., 0.]]),
    np.array([[0., 5.04892, 2.80194, 0., 0., -0.356896, 0., 6.49396, 0.], [0., 1., 
  0., 0., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0., 0., 
  0.], [0., 2.80194, 5.04892, 0., 0., -0.356896, 0., 6.49396, 
  0.], [0., 7.2959, 9.09783, 0., 0., -0.801938, 0., 14.5918, 0.], [0.,
   10.0978, 10.0978, 0., 0., -1., 0., 18.1957, 0.], [0., 9.09783, 
  7.2959, 0., 0., -0.801938, 0., 14.5918, 0.], [0., 0., 0., 0., 0., 
  0., 0., 1., 0.], [0., 2.34183, 2.34183, 0., 0., -0.219833, 0., 3., 
  0.]]),
    np.array([[0., 0., 7.2959, 10.5429, -2.24698, 0., 0., 14.5918, 0.], [0., 0., 
  4.24698, 4.24698, -1., 0., 0., 6.49396, 0.], [0., 0., 1., 0., 0., 
  0., 0., 0., 0.], [0., 0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 2., 
  6.49396, -1., 0., 0., 6.49396, 0.], [0., 0., 5.49396, 
  12.3448, -2.24698, 0., 0., 14.5918, 0.], [0., 0., 7.85086, 
  14.1468, -2.80194, 0., 0., 18.1957, 0.], [0., 0., 0., 0., 0., 0., 
  0., 1., 0.], [0., 0., 1.84787, 3.23191, -0.615957, 0., 0., 3., 0.]]),
    np.array([[0., 0., -2.80194, 14.1468, 7.85086, 0., 0., 18.1957, 0.], [0., 
  0., -2.24698, 12.3448, 5.49396, 0., 0., 14.5918, 0.], [0., 0., -1., 
  6.49396, 2., 0., 0., 6.49396, 0.], [0., 0., 0., 1., 0., 0., 0., 0., 
  0.], [0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., -1., 4.24698, 
  4.24698, 0., 0., 6.49396, 0.], [0., 0., -2.24698, 10.5429, 7.2959, 
  0., 0., 14.5918, 0.], [0., 0., 0., 0., 0., 0., 0., 1., 0.], [0., 
  0., -0.615957, 3.23191, 1.84787, 0., 0., 3., 0.]]),
    np.array([[0., -0.801938, 0., 0., 7.2959, 9.09783, 0., 14.5918, 0.], [0., -1., 
  0., 0., 10.0978, 10.0978, 0., 18.1957, 0.], [0., -0.801938, 0., 0., 
  9.09783, 7.2959, 0., 14.5918, 0.], [0., -0.356896, 0., 0., 5.04892, 
  2.80194, 0., 6.49396, 0.], [0., 0., 0., 0., 1., 0., 0., 0., 
  0.], [0., 0., 0., 0., 0., 1., 0., 0., 0.], [0., -0.356896, 0., 0., 
  2.80194, 5.04892, 0., 6.49396, 0.], [0., 0., 0., 0., 0., 0., 0., 1.,
   0.], [0., -0.219833, 0., 0., 2.34183, 2.34183, 0., 3., 0.]]),
    np.array([[0., 0., 0., -0.445042, 0., 3.24698, 4.69202, 6.49396, 0.], [0., 0., 
  0., -1., 0., 8.2959, 8.2959, 14.5918, 0.], [0., 0., 0., -1.24698, 
  0., 11.3448, 9.09783, 18.1957, 0.], [0., 0., 0., -1., 0., 10.0978, 
  6.49396, 14.5918, 0.], [0., 0., 0., -0.445042, 0., 5.49396, 2.44504,
   6.49396, 0.], [0., 0., 0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 
  0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 1., 0.], [0., 
  0., 0., -0.274127, 0., 2.61596, 2.122, 3., 0.]]),
    np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.], [4.24698, 0., 0., 0., 0., -1., 
  4.24698, 6.49396, 0.], [7.2959, 0., 0., 0., 0., -2.24698, 10.5429, 
  14.5918, 0.], [7.85086, 0., 0., 0., 0., -2.80194, 14.1468, 18.1957, 
  0.], [5.49396, 0., 0., 0., 0., -2.24698, 12.3448, 14.5918, 0.], [2.,
   0., 0., 0., 0., -1., 6.49396, 6.49396, 0.], [0., 0., 0., 0., 0., 
  0., 1., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 1., 0.], [1.84787, 0.,
   0., 0., 0., -0.615957, 3.23191, 3., 0.]]),
    np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0., 
  0., 0.], [2.80194, 5.04892, 0., 0., -0.356896, 0., 0., 0., 
  6.49396], [7.2959, 9.09783, 0., 0., -0.801938, 0., 0., 0., 
  14.5918], [10.0978, 10.0978, 0., 0., -1., 0., 0., 0., 
  18.1957], [9.09783, 7.2959, 0., 0., -0.801938, 0., 0., 0., 
  14.5918], [5.04892, 2.80194, 0., 0., -0.356896, 0., 0., 0., 
  6.49396], [2.34183, 2.34183, 0., 0., -0.219833, 0., 0., 0., 
  3.], [0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[0., 4.69202, 3.24698, 0., -0.445042, 0., 0., 0., 6.49396], [0., 1., 
  0., 0., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0., 0., 
  0.], [0., 2.44504, 5.49396, 0., -0.445042, 0., 0., 0., 
  6.49396], [0., 6.49396, 10.0978, 0., -1., 0., 0., 0., 14.5918], [0.,
   9.09783, 11.3448, 0., -1.24698, 0., 0., 0., 18.1957], [0., 8.2959, 
  8.2959, 0., -1., 0., 0., 0., 14.5918], [0., 2.122, 2.61596, 
  0., -0.274127, 0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 
  1.]]),
    np.array([[0., 0., 9.09783, 7.2959, 0., 0., -0.801938, 0., 14.5918], [0., 0., 
  5.04892, 2.80194, 0., 0., -0.356896, 0., 6.49396], [0., 0., 1., 0., 
  0., 0., 0., 0., 0.], [0., 0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 
  2.80194, 5.04892, 0., 0., -0.356896, 0., 6.49396], [0., 0., 7.2959, 
  9.09783, 0., 0., -0.801938, 0., 14.5918], [0., 0., 10.0978, 10.0978,
   0., 0., -1., 0., 18.1957], [0., 0., 2.34183, 2.34183, 0., 
  0., -0.219833, 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[0., -1.24698, 0., 11.3448, 9.09783, 0., 0., 0., 18.1957], [0., -1., 
  0., 10.0978, 6.49396, 0., 0., 0., 14.5918], [0., -0.445042, 0., 
  5.49396, 2.44504, 0., 0., 0., 6.49396], [0., 0., 0., 1., 0., 0., 0.,
   0., 0.], [0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., -0.445042, 0., 
  3.24698, 4.69202, 0., 0., 0., 6.49396], [0., -1., 0., 8.2959, 
  8.2959, 0., 0., 0., 14.5918], [0., -0.274127, 0., 2.61596, 2.122, 
  0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[-1., 0., 0., 0., 6.49396, 10.0978, 0., 0., 14.5918], [-1.24698, 0., 
  0., 0., 9.09783, 11.3448, 0., 0., 18.1957], [-1., 0., 0., 0., 
  8.2959, 8.2959, 0., 0., 14.5918], [-0.445042, 0., 0., 0., 4.69202, 
  3.24698, 0., 0., 6.49396], [0., 0., 0., 0., 1., 0., 0., 0., 
  0.], [0., 0., 0., 0., 0., 1., 0., 0., 0.], [-0.445042, 0., 0., 0., 
  2.44504, 5.49396, 0., 0., 6.49396], [-0.274127, 0., 0., 0., 2.122, 
  2.61596, 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[0., 0., -0.356896, 0., 0., 2.80194, 5.04892, 0., 6.49396], [0., 
  0., -0.801938, 0., 0., 7.2959, 9.09783, 0., 14.5918], [0., 0., -1., 
  0., 0., 10.0978, 10.0978, 0., 18.1957], [0., 0., -0.801938, 0., 0., 
  9.09783, 7.2959, 0., 14.5918], [0., 0., -0.356896, 0., 0., 5.04892, 
  2.80194, 0., 6.49396], [0., 0., 0., 0., 0., 1., 0., 0., 0.], [0., 
  0., 0., 0., 0., 0., 1., 0., 0.], [0., 0., -0.219833, 0., 0., 
  2.34183, 2.34183, 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
    np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0.], [4.24698, 0., 0., 0., 0., -1., 
  4.24698, 0., 6.49396], [7.2959, 0., 0., 0., 0., -2.24698, 10.5429, 
  0., 14.5918], [7.85086, 0., 0., 0., 0., -2.80194, 14.1468, 0., 
  18.1957], [5.49396, 0., 0., 0., 0., -2.24698, 12.3448, 0., 
  14.5918], [2., 0., 0., 0., 0., -1., 6.49396, 0., 6.49396], [0., 0., 
  0., 0., 0., 0., 1., 0., 0.], [1.84787, 0., 0., 0., 0., -0.615957, 
  3.23191, 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 1.]])]

    root = Node([3.30476487096249, 3.30476487096249, 3.30476487096249, \
3.30476487096249, 3.30476487096249, 3.30476487096249, \
3.30476487096249, 2.53284323061569, -1.00000000000000], [], words, False)
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

    m = 100
    print("maximum:",m)
    words = {}
    matrix = constructMatrix(words,m)
    #print(matrix)
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*13)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()