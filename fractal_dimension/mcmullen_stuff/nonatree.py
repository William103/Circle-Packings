from nonanode import Node
import math
import numpy as np
import time
from scipy.sparse import csr_matrix

def functions(n, z):
    duals = [[0.490290596565702+0.178451183290535j, 
  0.178451183290535], [0.260878177459588+0.451854257945975j, 
  0.178451183290535], [-0.0906020402178549+0.513829703507793j, 
  0.178451183290535], [-0.399688556347847+0.335378520217258j, 
  0.178451183290535], [-0.521756354919175, 
  0.178451183290535], [-0.399688556347847-0.335378520217258j, 
  0.178451183290535], [-0.0906020402178549-0.513829703507793j, 
  0.178451183290535], [0.260878177459588-0.451854257945975j, 
  0.178451183290535], [0.490290596565702-0.178451183290535j, 
  0.178451183290535], [1.00000000000000+0.363970234266202j, 
  0.363970234266202], [0.532088886237956+0.921604985106876j, 
  0.363970234266202], [-0.184792530904095+1.04801052091754j, 
  0.363970234266202], [-0.815207469095905+0.684040286651337j, 
  0.363970234266202], [-1.06417777247591, 
  0.363970234266202], [-0.815207469095905-0.684040286651337j, 
  0.363970234266202], [-0.184792530904095-1.04801052091754j, 
  0.363970234266202], [0.532088886237956-0.921604985106876j, 
  0.363970234266202], [1.00000000000000-0.363970234266202j, 
  0.363970234266202]]
    return duals[n][0] + np.conj((duals[n][1])**2 / (z-duals[n][0]))

def derivatives(n, z):
    duals = [[0.490290596565702 + 0.178451183290535j, 
  0.178451183290535], [0.260878177459588 + 0.451854257945975j, 
  0.178451183290535], [-0.0906020402178549 + 0.513829703507793j, 
  0.178451183290535], [-0.399688556347847 + 0.335378520217258j, 
  0.178451183290535], [-0.521756354919175, 
  0.178451183290535], [-0.399688556347847 - 0.335378520217258j, 
  0.178451183290535], [-0.0906020402178549 - 0.513829703507793j, 
  0.178451183290535], [0.260878177459588 - 0.451854257945975j, 
  0.178451183290535], [0.490290596565702 - 0.178451183290535j, 
  0.178451183290535], [1.00000000000000 + 0.363970234266202j, 
  0.363970234266202], [0.532088886237956 + 0.921604985106876j, 
  0.363970234266202], [-0.184792530904095 + 1.04801052091754j, 
  0.363970234266202], [-0.815207469095905 + 0.684040286651337j, 
  0.363970234266202], [-1.06417777247591, 
  0.363970234266202], [-0.815207469095905 - 0.684040286651337j, 
  0.363970234266202], [-0.184792530904095 - 1.04801052091754j, 
  0.363970234266202], [0.532088886237956 - 0.921604985106876j, 
  0.363970234266202], [1.00000000000000 - 0.363970234266202j, 
  0.363970234266202]]
    return abs(-1*((z-duals[n][0])**2) / ((duals[n][1])**2))

def samplePoint(word):
    points = [0.5148051263939873 + 0.18737374245506147j, 0.2739220863325669 + 
 0.47444697084327425j, -0.09513214222874763 + 
 0.5395211886831826j, -0.4196729841652397 + 
 0.3521474462281211j, -0.5478441726651339, -0.4196729841652397 - 
 0.3521474462281211j, -0.09513214222874763 - 
 0.5395211886831826j, 0.2739220863325669 - 
 0.47444697084327425j, 0.5148051263939873 - 
 0.18737374245506147j, 0.8 + 
 0.29117618741296186j, 0.42567110899036487 + 
 0.7372839880855011j, -0.1478340247232763 + 
 0.838408416734032j, -0.6521659752767237 + 
 0.54723222932107j, -0.8513422179807297, -0.6521659752767237 - 
 0.54723222932107j, -0.1478340247232763 - 
 0.838408416734032j, 0.42567110899036487 - 
 0.7372839880855011j, 0.8 - 0.29117618741296186j]
    p = points[word[-1]]
    for letter in word[-2::-1]:
        p = functions(letter, p)
    return p

def sampleValue(word):
    return derivatives(word[0], samplePoint(word))

def generateTree(words, dc):
    generators = [np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 
  0., 0., 0., 0., 0., 0.], [2.87939, 5.41147, 0., 0., 0., -0.226682, 
  0., 0., 0., 7.06418, 0.], [8.29086, 11.1702, 0., 0., 0., -0.573978, 
  0., 0., 0., 17.8871, 0.], [13.7023, 15.5817, 0., 0., 0., -0.879385, 
  0., 0., 0., 27.4047, 0.], [16.5817, 16.5817, 0., 0., 0., -1., 0., 
  0., 0., 31.1634, 0.], [15.5817, 13.7023, 0., 0., 0., -0.879385, 0., 
  0., 0., 27.4047, 0.], [11.1702, 8.29086, 0., 0., 0., -0.573978, 0., 
  0., 0., 17.8871, 0.], [5.41147, 2.87939, 0., 0., 0., -0.226682, 0., 
  0., 0., 7.06418, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 
  0.], [2.19665, 2.19665, 0., 0., 0., -0.128356, 0., 0., 0., 3., 0.]]), 
  np.array([[0., 4.92702, 3.53209, 0., -0.394931, 0., 0., 0., 0., 7.06418, 
  0.], [0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 
  0., 0., 0., 0., 0., 0., 0.], [0., 2.39493, 6.06418, 0., -0.394931, 
  0., 0., 0., 0., 7.06418, 0.], [0., 7.06418, 12.8229, 0., -1., 0., 
  0., 0., 0., 17.8871, 0.], [0., 11.8229, 18.1138, 0., -1.53209, 0., 
  0., 0., 0., 27.4047, 0.], [0., 14.4446, 19.4611, 0., -1.74223, 0., 
  0., 0., 0., 31.1634, 0.], [0., 13.7023, 16.2344, 0., -1.53209, 0., 
  0., 0., 0., 27.4047, 0.], [0., 9.94356, 9.94356, 0., -1., 0., 0., 
  0., 0., 17.8871, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 
  0.], [0., 1.92234, 2.56624, 0., -0.223625, 0., 0., 0., 0., 3., 0.]]),
  np.array([[0., 0., 11.1702, 8.29086, 0., 0., 0., -0.573978, 0., 17.8871, 
  0.], [0., 0., 5.41147, 2.87939, 0., 0., 0., -0.226682, 0., 7.06418, 
  0.], [0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 1., 
  0., 0., 0., 0., 0., 0., 0.], [0., 0., 2.87939, 5.41147, 0., 0., 
  0., -0.226682, 0., 7.06418, 0.], [0., 0., 8.29086, 11.1702, 0., 0., 
  0., -0.573978, 0., 17.8871, 0.], [0., 0., 13.7023, 15.5817, 0., 0., 
  0., -0.879385, 0., 27.4047, 0.], [0., 0., 16.5817, 16.5817, 0., 0., 
  0., -1., 0., 31.1634, 0.], [0., 0., 15.5817, 13.7023, 0., 0., 
  0., -0.879385, 0., 27.4047, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 
  0., 1., 0.], [0., 0., 2.19665, 2.19665, 0., 0., 0., -0.128356, 0., 
  3., 0.]]),
  np.array([[0., -1.53209, 0., 18.1138, 11.8229, 0., 0., 0., 0., 27.4047, 
  0.], [0., -1., 0., 12.8229, 7.06418, 0., 0., 0., 0., 17.8871, 
  0.], [0., -0.394931, 0., 6.06418, 2.39493, 0., 0., 0., 0., 7.06418, 
  0.], [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 
  1., 0., 0., 0., 0., 0., 0.], [0., -0.394931, 0., 3.53209, 4.92702, 
  0., 0., 0., 0., 7.06418, 0.], [0., -1., 0., 9.94356, 9.94356, 0., 
  0., 0., 0., 17.8871, 0.], [0., -1.53209, 0., 16.2344, 13.7023, 0., 
  0., 0., 0., 27.4047, 0.], [0., -1.74223, 0., 19.4611, 14.4446, 0., 
  0., 0., 0., 31.1634, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 
  0.], [0., -0.223625, 0., 2.56624, 1.92234, 0., 0., 0., 0., 3., 0.]]),
  np.array([[0., -1.13716, 0., 0., 17.7189, 15.5817, 0., 0., 0., 31.1634, 
  0.], [0., -1., 0., 0., 16.5817, 12.8229, 0., 0., 0., 27.4047, 
  0.], [0., -0.652704, 0., 0., 11.8229, 7.71688, 0., 0., 0., 17.8871, 
  0.], [0., -0.257773, 0., 0., 5.66925, 2.6527, 0., 0., 0., 7.06418, 
  0.], [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 
  0., 1., 0., 0., 0., 0., 0.], [0., -0.257773, 0., 0., 3.13716, 
  5.18479, 0., 0., 0., 7.06418, 0.], [0., -0.652704, 0., 0., 8.94356, 
  10.5963, 0., 0., 0., 17.8871, 0.], [0., -1., 0., 0., 14.7023, 
  14.7023, 0., 0., 0., 27.4047, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 
  0., 1., 0.], [0., -0.145961, 0., 0., 2.34261, 2.0683, 0., 0., 0., 
  3., 0.]]),
  np.array([[0., 0., 0., 0., 0., 10.2909, 21.9932, -3.87939, 0., 27.4047, 
  0.], [0., 0., 0., 0., 0., 12.7023, 23.8726, -4.41147, 0., 31.1634, 
  0.], [0., 0., 0., 0., 0., 12.1702, 20.1138, -3.87939, 0., 27.4047, 
  0.], [0., 0., 0., 0., 0., 8.94356, 12.4757, -2.53209, 0., 17.8871, 
  0.], [0., 0., 0., 0., 0., 4.53209, 4.53209, -1., 0., 7.06418, 
  0.], [0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 
  0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 2., 7.06418, -1., 
  0., 7.06418, 0.], [0., 0., 0., 0., 0., 6.06418, 15.355, -2.53209, 
  0., 17.8871, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.], [0.,
   0., 0., 0., 0., 1.69871, 3.13247, -0.566237, 0., 3., 0.]]),
   np.array([[-1., 0., 0., 0., 0., 0., 7.06418, 12.8229, 0., 17.8871, 
  0.], [-1.53209, 0., 0., 0., 0., 0., 11.8229, 18.1138, 0., 27.4047, 
  0.], [-1.74223, 0., 0., 0., 0., 0., 14.4446, 19.4611, 0., 31.1634, 
  0.], [-1.53209, 0., 0., 0., 0., 0., 13.7023, 16.2344, 0., 27.4047, 
  0.], [-1., 0., 0., 0., 0., 0., 9.94356, 9.94356, 0., 17.8871, 
  0.], [-0.394931, 0., 0., 0., 0., 0., 4.92702, 3.53209, 0., 7.06418, 
  0.], [0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 0., 
  0., 0., 0., 1., 0., 0., 0.], [-0.394931, 0., 0., 0., 0., 0., 
  2.39493, 6.06418, 0., 7.06418, 0.], [0., 0., 0., 0., 0., 0., 0., 0.,
   0., 1., 0.], [-0.223625, 0., 0., 0., 0., 0., 1.92234, 2.56624, 0., 
  3., 0.]]),
  np.array([[0., 0., 0., 0., 0., 0., 0., 3., 5.53209, 5.29813, -1.76604], [0., 
  0., 0., 0., 0., 0., 0., 8.59627, 11.4757, 13.4153, -4.47178], [0., 
  0., 0., 0., 0., 0., 0., 14.1702, 16.0496, 20.5535, -6.85117], [0., 
  0., 0., 0., 0., 0., 0., 17.1138, 17.1138, 23.3726, -7.79086], [0., 
  0., 0., 0., 0., 0., 0., 16.0496, 14.1702, 20.5535, -6.85117], [0., 
  0., 0., 0., 0., 0., 0., 11.4757, 8.59627, 13.4153, -4.47178], [0., 
  0., 0., 0., 0., 0., 0., 5.53209, 3., 5.29813, -1.76604], [0., 0., 
  0., 0., 0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 
  0., 1., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.], [0., 
  0., 0., 0., 0., 0., 0., 2.26495, 2.26495, 2., -1.]]),
  np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [7.06418, -1., 0., 0., 
  0., 0., 0., 0., 2., 7.06418, 0.], [15.355, -2.53209, 0., 0., 0., 0.,
   0., 0., 6.06418, 17.8871, 0.], [21.9932, -3.87939, 0., 0., 0., 0., 
  0., 0., 10.2909, 27.4047, 0.], [23.8726, -4.41147, 0., 0., 0., 0., 
  0., 0., 12.7023, 31.1634, 0.], [20.1138, -3.87939, 0., 0., 0., 0., 
  0., 0., 12.1702, 27.4047, 0.], [12.4757, -2.53209, 0., 0., 0., 0., 
  0., 0., 8.94356, 17.8871, 0.], [4.53209, -1., 0., 0., 0., 0., 0., 
  0., 4.53209, 7.06418, 0.], [0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 
  0.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 
  0.], [3.13247, -0.566237, 0., 0., 0., 0., 0., 0., 1.69871, 3., 0.]]),
  np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 
  0., 0., 0., 0., 0., 0.], [2.87939, 5.41147, 0., 0., 0., -0.226682, 
  0., 0., 0., 0., 7.06418], [8.29086, 11.1702, 0., 0., 0., -0.573978, 
  0., 0., 0., 0., 17.8871], [13.7023, 15.5817, 0., 0., 0., -0.879385, 
  0., 0., 0., 0., 27.4047], [16.5817, 16.5817, 0., 0., 0., -1., 0., 
  0., 0., 0., 31.1634], [15.5817, 13.7023, 0., 0., 0., -0.879385, 0., 
  0., 0., 0., 27.4047], [11.1702, 8.29086, 0., 0., 0., -0.573978, 0., 
  0., 0., 0., 17.8871], [5.41147, 2.87939, 0., 0., 0., -0.226682, 0., 
  0., 0., 0., 7.06418], [2.19665, 2.19665, 0., 0., 0., -0.128356, 0., 
  0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
  np.array([[0., 4.92702, 3.53209, 0., -0.394931, 0., 0., 0., 0., 0., 
  7.06418], [0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 1.,
   0., 0., 0., 0., 0., 0., 0., 0.], [0., 2.39493, 6.06418, 
  0., -0.394931, 0., 0., 0., 0., 0., 7.06418], [0., 7.06418, 12.8229, 
  0., -1., 0., 0., 0., 0., 0., 17.8871], [0., 11.8229, 18.1138, 
  0., -1.53209, 0., 0., 0., 0., 0., 27.4047], [0., 14.4446, 19.4611, 
  0., -1.74223, 0., 0., 0., 0., 0., 31.1634], [0., 13.7023, 16.2344, 
  0., -1.53209, 0., 0., 0., 0., 0., 27.4047], [0., 9.94356, 9.94356, 
  0., -1., 0., 0., 0., 0., 0., 17.8871], [0., 1.92234, 2.56624, 
  0., -0.223625, 0., 0., 0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0.,
   0., 0., 0., 1.]]),
   np.array([[0., 0., 8.94356, 12.4757, -2.53209, 0., 0., 0., 0., 0., 
  17.8871], [0., 0., 4.53209, 4.53209, -1., 0., 0., 0., 0., 0., 
  7.06418], [0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 0.,
   1., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 2., 7.06418, -1., 0., 0.,
   0., 0., 0., 7.06418], [0., 0., 6.06418, 15.355, -2.53209, 0., 0., 
  0., 0., 0., 17.8871], [0., 0., 10.2909, 21.9932, -3.87939, 0., 0., 
  0., 0., 0., 27.4047], [0., 0., 12.7023, 23.8726, -4.41147, 0., 0., 
  0., 0., 0., 31.1634], [0., 0., 12.1702, 20.1138, -3.87939, 0., 0., 
  0., 0., 0., 27.4047], [0., 0., 1.69871, 3.13247, -0.566237, 0., 0., 
  0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
  np.array([[0., 0., -3.87939, 21.9932, 10.2909, 0., 0., 0., 0., 0., 
  27.4047], [0., 0., -2.53209, 15.355, 6.06418, 0., 0., 0., 0., 0., 
  17.8871], [0., 0., -1., 7.06418, 2., 0., 0., 0., 0., 0., 
  7.06418], [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 0.,
   0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., -1., 4.53209, 4.53209, 
  0., 0., 0., 0., 0., 7.06418], [0., 0., -2.53209, 12.4757, 8.94356, 
  0., 0., 0., 0., 0., 17.8871], [0., 0., -3.87939, 20.1138, 12.1702, 
  0., 0., 0., 0., 0., 27.4047], [0., 0., -4.41147, 23.8726, 12.7023, 
  0., 0., 0., 0., 0., 31.1634], [0., 0., -0.566237, 3.13247, 1.69871, 
  0., 0., 0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
  1.]]),
  np.array([[0., 0., 0., 0., 12.7023, 23.8726, -4.41147, 0., 0., 0., 
  31.1634], [0., 0., 0., 0., 12.1702, 20.1138, -3.87939, 0., 0., 0., 
  27.4047], [0., 0., 0., 0., 8.94356, 12.4757, -2.53209, 0., 0., 0., 
  17.8871], [0., 0., 0., 0., 4.53209, 4.53209, -1., 0., 0., 0., 
  7.06418], [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., 0.,
   0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 2., 7.06418, -1.,
   0., 0., 0., 7.06418], [0., 0., 0., 0., 6.06418, 15.355, -2.53209, 
  0., 0., 0., 17.8871], [0., 0., 0., 0., 10.2909, 21.9932, -3.87939, 
  0., 0., 0., 27.4047], [0., 0., 0., 0., 1.69871, 3.13247, -0.566237, 
  0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
  np.array([[0., 0., -1., 0., 0., 14.7023, 14.7023, 0., 0., 0., 27.4047], [0., 
  0., -1.13716, 0., 0., 17.7189, 15.5817, 0., 0., 0., 31.1634], [0., 
  0., -1., 0., 0., 16.5817, 12.8229, 0., 0., 0., 27.4047], [0., 
  0., -0.652704, 0., 0., 11.8229, 7.71688, 0., 0., 0., 17.8871], [0., 
  0., -0.257773, 0., 0., 5.66925, 2.6527, 0., 0., 0., 7.06418], [0., 
  0., 0., 0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0., 
  1., 0., 0., 0., 0.], [0., 0., -0.257773, 0., 0., 3.13716, 5.18479, 
  0., 0., 0., 7.06418], [0., 0., -0.652704, 0., 0., 8.94356, 10.5963, 
  0., 0., 0., 17.8871], [0., 0., -0.145961, 0., 0., 2.34261, 2.0683, 
  0., 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]]),
  np.array([[-1., 0., 0., 0., 0., 0., 7.06418, 12.8229, 0., 0., 
  17.8871], [-1.53209, 0., 0., 0., 0., 0., 11.8229, 18.1138, 0., 0., 
  27.4047], [-1.74223, 0., 0., 0., 0., 0., 14.4446, 19.4611, 0., 0., 
  31.1634], [-1.53209, 0., 0., 0., 0., 0., 13.7023, 16.2344, 0., 0., 
  27.4047], [-1., 0., 0., 0., 0., 0., 9.94356, 9.94356, 0., 0., 
  17.8871], [-0.394931, 0., 0., 0., 0., 0., 4.92702, 3.53209, 0., 0., 
  7.06418], [0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., 0.,
   0., 0., 0., 0., 1., 0., 0., 0.], [-0.394931, 0., 0., 0., 0., 0., 
  2.39493, 6.06418, 0., 0., 7.06418], [-0.223625, 0., 0., 0., 0., 0., 
  1.92234, 2.56624, 0., 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 
  0., 1.]]),
  np.array([[0., 0., 0., -0.226682, 0., 0., 0., 2.87939, 5.41147, 0., 
  7.06418], [0., 0., 0., -0.573978, 0., 0., 0., 8.29086, 11.1702, 0., 
  17.8871], [0., 0., 0., -0.879385, 0., 0., 0., 13.7023, 15.5817, 0., 
  27.4047], [0., 0., 0., -1., 0., 0., 0., 16.5817, 16.5817, 0., 
  31.1634], [0., 0., 0., -0.879385, 0., 0., 0., 15.5817, 13.7023, 0., 
  27.4047], [0., 0., 0., -0.573978, 0., 0., 0., 11.1702, 8.29086, 0., 
  17.8871], [0., 0., 0., -0.226682, 0., 0., 0., 5.41147, 2.87939, 0., 
  7.06418], [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.], [0., 0., 0.,
   0., 0., 0., 0., 0., 1., 0., 0.], [0., 0., 0., -0.128356, 0., 0., 
  0., 2.19665, 2.19665, 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 
  0., 1.]]),
  np.array([[1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], [5.18479, 0., 0., 0., 
  0., -0.257773, 0., 0., 3.13716, 0., 7.06418], [10.5963, 0., 0., 0., 
  0., -0.652704, 0., 0., 8.94356, 0., 17.8871], [14.7023, 0., 0., 0., 
  0., -1., 0., 0., 14.7023, 0., 27.4047], [15.5817, 0., 0., 0., 
  0., -1.13716, 0., 0., 17.7189, 0., 31.1634], [12.8229, 0., 0., 0., 
  0., -1., 0., 0., 16.5817, 0., 27.4047], [7.71688, 0., 0., 0., 
  0., -0.652704, 0., 0., 11.8229, 0., 17.8871], [2.6527, 0., 0., 0., 
  0., -0.257773, 0., 0., 5.66925, 0., 7.06418], [0., 0., 0., 0., 0., 
  0., 0., 0., 1., 0., 0.], [2.0683, 0., 0., 0., 0., -0.145961, 0., 0.,
   2.34261, 0., 3.], [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.]])]

    root = Node([3.92380440016309, 3.92380440016309, 3.92380440016309, \
3.92380440016309, 3.92380440016309, 3.92380440016309, \
3.92380440016309, 3.92380440016309, 3.92380440016309, \
2.03960672916147, -1.00000000000000], [], words, False)
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
    m = 100
    print("maximum:",m)
    matrix = constructMatrix(words, m)
    #print(matrix)
    
    print("construction (s): %f\nconstruction (m): %f" % (time.time()-start, (time.time()-start)/60))
    #print(csr_matrix.transpose(matrix))
    print(csr_matrix.count_nonzero(matrix), (matrix.shape[0])*15)
    secantMethod(matrix,matrix.shape[0],1,1.33,1.34,10**(-10),1000)

    print("total (s): %f\ntotal (m): %f" % (time.time()-start, (time.time()-start)/60))

main()