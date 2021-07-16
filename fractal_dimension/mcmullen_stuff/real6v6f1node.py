import numpy as np
import math

class Node:
    def __init__(self, tuple, word, words, infertile):
        self.tuple = tuple
        self.children = []
        self.word = word
        self.infertile = infertile
        words[str(word)] = self

    def dc_not_too_big(self, g, generator, dual_curvature):
        temp = np.matmul(generator, self.tuple)
        if g == 0:
            temp = [temp[i] for i in (2,3,4,5)]
            #print(temp)
            if temp[-1]<0:
                return True
            if temp[0]*temp[1]*temp[2] + temp[0]*temp[1]*temp[3] + temp[0]*temp[2]*temp[3] + temp[1]*temp[2]*temp[3] < 0:
                return False
            else:
                return math.sqrt(temp[0]*temp[1]*temp[2] + temp[0]*temp[1]*temp[3] + temp[0]*temp[2]*temp[3] + temp[1]*temp[2]*temp[3])/math.sqrt(temp[0] + temp[1] + temp[2] + temp[3]) < dual_curvature
        elif g == 1:
            temp = [temp[i] for i in (1,2,5)]
        elif g == 2:
            temp = [temp[i] for i in (0,1,5)]
        elif g == 3:
            temp = [temp[i] for i in (0,4,5)]
        elif g == 4:
            temp = [temp[i] for i in (0,1,2,3)]
            #print(temp)
            if temp[0]*temp[1]*temp[2] + temp[0]*temp[1]*temp[3] + temp[0]*temp[2]*temp[3] + temp[1]*temp[2]*temp[3] < 0:
                return False
            return math.sqrt(temp[0]*temp[1]*temp[2] + temp[0]*temp[1]*temp[3] + temp[0]*temp[2]*temp[3] + temp[1]*temp[2]*temp[3])/math.sqrt(temp[0] + temp[1] + temp[2] + temp[3]) < dual_curvature
        elif g == 5:
            temp = [temp[i] for i in (0,3,4)]
        if temp[0]*temp[1]+temp[0]*temp[2]+temp[1]*temp[2] < 0:
            return False
        return math.sqrt(temp[0]*temp[1]+temp[0]*temp[2]+temp[1]*temp[2])<dual_curvature
    
    def next_generation(self, words, dual_curvature, generators):
        if self.infertile:
            return [self]
        if self.children == []:
            for g,generator in enumerate(generators):
                if len(self.word) == 0 or g != self.word[-1]:
                    self.children.append(Node(np.matmul(generator,self.tuple), self.word[:] + [g], words, not self.dc_not_too_big(g, generator, dual_curvature)))
        if self.children == []:
            return [self]
        else:
            return self.children
    
    def next(self):
        if self.children == []:
            return [self]
        else:
            return self.children
    
    def leaves(self):
        current_leaves = [self]
        while True:
            new_leaves = []
            for leaf in current_leaves:
                new_leaves += leaf.next()
            if current_leaves == new_leaves:
                break
            else:
                current_leaves = new_leaves
        return current_leaves