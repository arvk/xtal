import numpy as np

class AtTraj:
    # Define variables
    box = np.ndarray((3,3))
    basisa = box[0,:]
    basisb = box[1,:]
    basisc = box[2,:]

    def __init__(self):
        print 'Atomic trajectory initialized'
