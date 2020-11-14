from xtal import xtal
import numpy as np

u = xtal.AtTraj()
u.read_trajectory_vasp('XDATCAR')
u.unwrap_trajectory()

dt = 10.0 # timestep between frames in a trajectory
corr_length = 6000  # 6000 fs
offset = 300        # 300 fs

print(len(u.snaplist))

atomIDs = np.arange(len(u.snaplist[0].atomlist))
num_atoms = len(u.snaplist[0].atomlist)

timeIDs = np.arange(0, np.floor(offset/dt),dtype=np.int)
actual_time = timeIDs * dt
msd_array = np.zeros_like(timeIDs, dtype=np.float)

num_initial = len(np.arange(0, len(u.snaplist)-np.ceil(corr_length/dt) , np.ceil(offset/dt), dtype=np.int))

for initframe in np.arange(0, len(u.snaplist)-np.ceil(corr_length/dt) , np.ceil(offset/dt), dtype=np.int):
    for delta in np.arange(0, np.floor(offset/dt),dtype=np.int):
        msd = 0.0
        for ID in atomIDs:
             dist = u.snaplist[initframe+delta].atomlist[ID].cart - u.snaplist[initframe].atomlist[ID].cart
             msd += (np.linalg.norm(dist)*np.linalg.norm(dist))
        # print(msd)
        # print(num_atoms)
        # print(num_initial)
        # print(msd / (num_initial*num_atoms))
        msd_array[delta] += msd / (num_initial*num_atoms)

print(msd_array)

np.savetxt('msd.txt',np.vstack((actual_time,msd_array)).T, fmt='%9.6f')
