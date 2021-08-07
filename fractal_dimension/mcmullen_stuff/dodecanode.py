import numpy as np
import math
from random import sample

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
            temp = [temp[i] for i in (0,1,16,12,17)]
        elif g == 1:
            temp = [temp[i] for i in (0,1,13,12,14)]
        elif g == 2:
            temp = [temp[i] for i in (9,11,19,8,15)]
        elif g == 3:
            temp = [temp[i] for i in (3,7,19,17,12)]
        elif g == 4:
            temp = [temp[i] for i in (2,13,4,8,15)]
        elif g == 5:
            temp = [temp[i] for i in (4,13,1,16,18)]
        elif g == 6:
            temp = [temp[i] for i in (7,11,9,5,19)]
        elif g == 7:
            temp = [temp[i] for i in (3,7,11,10,6)]
        elif g == 8:
            temp = [temp[i] for i in (0,14,5,19,17)]
        elif g == 9:
            temp = [temp[i] for i in (3,6,18,16,12)]
        elif g == 10:
            temp = [temp[i] for i in (2,14,5,9,15)]
        elif g == 11:
            temp = [temp[i] for i in (4,8,10,6,18)]
        if temp[0] < 0:
            return math.sqrt(-((temp[0]*temp[1]*temp[2] + temp[0]*temp[1]*temp[3] + temp[0]*temp[2]*temp[3] - temp[1]*temp[2]*temp[3] + temp[0]*temp[1]*temp[4] + 
  temp[0]*temp[2]*temp[4] - temp[1]*temp[2]*temp[4] + temp[0]*temp[3]*temp[4] - temp[1]*temp[3]*temp[4] - temp[2]*temp[3]*temp[4] - 
  math.sqrt(-4*temp[0]*temp[1]*temp[2]*temp[3]*(temp[0] - temp[1] - temp[2] - temp[3] - temp[4])*temp[4] + (temp[2]*temp[3]*temp[4] + 
     temp[1]*(temp[3]*temp[4] + temp[2]*(temp[3] + temp[4])) - 
     temp[0]*(temp[3]*temp[4] + temp[2]*(temp[3] + temp[4]) + temp[1]*(temp[2] + temp[3] + temp[4])))**2))/(-temp[0] + temp[1] + 
  temp[2] + temp[3] + temp[4])))/math.sqrt(2) < dual_curvature
        else:
            return (1/math.sqrt(2))*(math.sqrt(temp[0]*temp[1] + temp[0]*temp[2] + temp[1]*temp[2] + temp[0]*temp[3] + temp[1]*temp[3] + temp[2]*temp[3] + 
    temp[0]*temp[4] + temp[1]*temp[4] + temp[2]*temp[4] + temp[3]*temp[4] - 
    math.sqrt(-4*temp[0]*temp[1]*temp[2]*temp[3] - 4*temp[0]*temp[1]*temp[2]*temp[4] - 4*temp[0]*temp[1]*temp[3]*temp[4] - 
     4*temp[0]*temp[2]*temp[3]*temp[4] - 
     4*temp[1]*temp[2]*temp[3]*temp[4] + (temp[2]*temp[3] + temp[2]*temp[4] + temp[3]*temp[4] + temp[1]*(temp[2] + temp[3] + temp[4]) + 
       temp[0]*(temp[1] + temp[2] + temp[3] + temp[4]))**2))) < dual_curvature

    def next_generation(self, words, dual_curvature, generators):
        if self.infertile:
            return [self]
        if self.children == []:
            for g,generator in enumerate(generators):
                if len(self.word) == 0 or g != self.word[-1]:
                    if self.word + [g] == [1, 8, 15, 8, 0] or self.word + [g] == [8, 15, 8, 0]:
                        print(self.tuple)
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