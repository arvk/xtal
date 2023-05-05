import xtal
from utils import make_polycrystal
import numpy as np

u = xtal.AtTraj()
u.read_snapshot_vasp('SINGLE_CRYSTAL_POSCAR')
v = make_polycrystal(u.snaplist[0], 20, np.array([100,100,100]), remove_duplicates = 1.0)
v.snaplist[0].write_snapshot_vasp('POLYCRYSTAL_POSCAR',False)
