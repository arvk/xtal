from periodictable import core, mass, density
import copy
import numpy as np

chem_ref = core.PeriodicTable(table="X")
mass.init(chem_ref)
density.init(chem_ref)
chem_ref.X = copy.deepcopy(chem_ref.n)
chem_ref.X.symbol = 'X'
chem_ref.Xa = copy.deepcopy(chem_ref.n)
chem_ref.Xa.symbol = 'Xa'
chem_ref.Xb = copy.deepcopy(chem_ref.n)
chem_ref.Xb.symbol = 'Xb'



def is_sierpinski_carpet_filled(level, coords):
    '''Calculate if the given fractional coordinates correspond to a filed pixel
    in the Sierpinksi carpet fractal of a given level

    Multiply the fractional coordinate with 3^n (for n = level..1) and check to see
    if it leaves a reminder of 1 upon division by 3 (i.e. it is the middle cell at any level)'''

    multiplier = 3**level

    x = int(coords[0]*multiplier) # pylint: disable=invalid-name
    y = int(coords[1]*multiplier) # pylint: disable=invalid-name

    while True:
        if x == 0 and y == 0:
            break

        if x%3 == 1 and y%3 == 1:
            return False

        x = int(x/3) # pylint: disable=invalid-name
        y = int(y/3) # pylint: disable=invalid-name

    return True


def vector_to_quaternion(vector):
    return np.array([0.0,vector[0],vector[1],vector[2]])

def quaternion_to_vector(q):
    return np.array(q[1:])

def q_multiply(q1,q2):
    '''Quaternion multiplication. Returns product of two input quaternions'''
    t0 = (q1[0]*q2[0]) - (q1[1]*q2[1]) - (q1[2]*q2[2]) - (q1[3]*q2[3])
    t1 = (q1[0]*q2[1]) + (q1[1]*q2[0]) - (q1[2]*q2[3]) + (q1[3]*q2[2])
    t2 = (q1[0]*q2[2]) + (q1[1]*q2[3]) + (q1[2]*q2[0]) - (q1[3]*q2[1])
    t3 = (q1[0]*q2[3]) - (q1[1]*q2[2]) + (q1[2]*q2[1]) + (q1[3]*q2[0])
    return np.array([t0,t1,t2,t3])

def q_inv(q):
    '''Quaternion inverse. Returns inverse of an input quaternion'''
    return np.array([q[0],-q[1],-q[2],-q[3]])
