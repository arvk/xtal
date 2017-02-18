#! /usr/bin/python
### Make Origami vertices using the Miura-Ori design

import numpy as np
from scipy.interpolate import griddata


def make_ori_MO(m_in,n_in):

# Origami parameters
    ori_size = 100.0
    m = m_in*2
    n = n_in*2
    Tx = 2.0*ori_size/(m)
    Ty = 2.0*ori_size/(n)
    hx = 1.0
    hy = 1.0

    Vx = np.zeros([m,n])
    Vy = np.zeros([m,n])
    Vz = np.zeros([m,n])

    for i in range(1,m+1):
        for j in range(1,n+1):
            Vx[i-1,j-1] = (i-1)*Tx/2
            Vy[i-1,j-1] = ((j-1)*Ty/2) + ((1 + (-1)**i)/2)*(hx/hy)*np.sqrt((Ty/2)**2 + hy**2)
            Vz[i-1,j-1] = (1 + (-1)**j)*hy/2




    # Make creases
    Vtx = np.zeros([m,n])
    Vty = np.zeros([m,n])

    Vtx[0,0] = 0.0
    Vty[0,0] = 0.0

    for i in range(1,m):
        vec_i_1 = np.array([ Vx[i-1,0], Vy[i-1,0], Vz[i-1,0] ])
        vec_i_2 = np.array([ Vx[i-1,1], Vy[i-1,2], Vz[i-1,2] ])
        vec_ip1_1 = np.array([ Vx[i,0], Vy[i,0], Vz[i,0] ])
        vec_ip1_2 = np.array([ Vx[i,1], Vy[i,1], Vz[i,1] ])

        cos_phi = ((np.linalg.norm(vec_i_1-vec_i_2)**2) + (np.linalg.norm(vec_i_1-vec_ip1_1)**2) - (np.linalg.norm(vec_i_2-vec_ip1_1)**2))/(2.0*np.linalg.norm(vec_i_1-vec_i_2)*np.linalg.norm(vec_i_1-vec_ip1_1))
        cos_chi = ((np.linalg.norm(vec_ip1_1-vec_ip1_2)**2) + (np.linalg.norm(vec_i_1-vec_ip1_1)**2) - (np.linalg.norm(vec_ip1_2-vec_i_1)**2))/(2.0*np.linalg.norm(vec_ip1_1-vec_ip1_2)*np.linalg.norm(vec_i_1-vec_ip1_1))


        if ((cos_phi <= 1.0) and (cos_phi >= 0.0)):
            phi = np.arccos(cos_phi)
            Vtx[i,0] = Vtx[i-1,0] + (np.sin(phi)*np.linalg.norm(vec_i_1-vec_ip1_1))
            Vty[i,0] = Vty[i-1,0] + (np.cos(phi)*np.linalg.norm(vec_i_1-vec_ip1_1))

        else:
            chi = np.arccos(cos_chi)
            Vtx[i,0] = Vtx[i-1,0] + (np.sin(chi)*np.linalg.norm(vec_i_1-vec_ip1_1))
            Vty[i,0] = Vty[i-1,0] - (np.cos(chi)*np.linalg.norm(vec_i_1-vec_ip1_1))



    for j in range(1,n):
        for i in range(0,m):
            vec_i_j = np.array([ Vx[i,j], Vy[i,j], Vz[i,j] ])
            vec_i_jm1 = np.array([ Vx[i,j-1], Vy[i,j-1], Vz[i,j-1] ])
            Vtx[i,j] = Vtx[i,j-1]
            Vty[i,j] = Vty[i,j-1] + np.linalg.norm(vec_i_jm1-vec_i_j)



    # Convert to tuples

    #initpoints = np.array([])
    #initpoints = np.append(Vtx.reshape(m*n,1),Vty.reshape(m*n,1),1)

    crease_x = Vtx.reshape(m*n,1)
    crease_y = Vty.reshape(m*n,1)

    initpoints_i_j = np.array([])
    initpoints_i_j = np.append(crease_x,crease_y,1)  # i,j
    initpoints_ip1_j = np.append(crease_x+ori_size,crease_y,1)  # i+1,j
    initpoints_im1_j = np.append(crease_x-ori_size,crease_y,1)  # i-1,j
    initpoints_i_jp1 = np.append(crease_x,crease_y+ori_size,1)  # i,j+1
    initpoints_i_jm1 = np.append(crease_x,crease_y-ori_size,1)  # i,j-1
    initpoints_ip1_jp1 = np.append(crease_x+ori_size,crease_y+ori_size,1)  # i+1,j+1
    initpoints_ip1_jm1 = np.append(crease_x+ori_size,crease_y-ori_size,1)  # i+1,j-1
    initpoints_im1_jp1 = np.append(crease_x-ori_size,crease_y+ori_size,1)  # i-1,j+1
    initpoints_im1_jm1 = np.append(crease_x-ori_size,crease_y-ori_size,1)  # i-1,j-1

    initpoints = initpoints_i_j
    initpoints = np.append(initpoints,initpoints_ip1_j,0)
    initpoints = np.append(initpoints,initpoints_im1_j,0)
    initpoints = np.append(initpoints,initpoints_i_jp1,0)
    initpoints = np.append(initpoints,initpoints_i_jm1,0)
    initpoints = np.append(initpoints,initpoints_ip1_jp1,0)
    initpoints = np.append(initpoints,initpoints_ip1_jm1,0)
    initpoints = np.append(initpoints,initpoints_im1_jp1,0)
    initpoints = np.append(initpoints,initpoints_im1_jm1,0)

    ## Normalize to [0,1]
    initpoints = initpoints/ori_size




    # Convert Origami vertices to tuples

    values_ij = Vx.reshape(m*n,1)
    values = values_ij
    values = np.append(values,values_ij+ori_size,0)
    values = np.append(values,values_ij-ori_size,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij+ori_size,0)
    values = np.append(values,values_ij+ori_size,0)
    values = np.append(values,values_ij-ori_size,0)
    values = np.append(values,values_ij-ori_size,0)
    valuesx = values/ori_size

    values_ij = Vy.reshape(m*n,1)
    values = values_ij
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij+ori_size,0)
    values = np.append(values,values_ij-ori_size,0)
    values = np.append(values,values_ij+ori_size,0)
    values = np.append(values,values_ij-ori_size,0)
    values = np.append(values,values_ij+ori_size,0)
    values = np.append(values,values_ij-ori_size,0)
    valuesy = values/ori_size

    values_ij = Vz.reshape(m*n,1)
    values = values_ij
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    values = np.append(values,values_ij,0)
    valuesz = values/ori_size


    folded = valuesx
    folded = np.append(folded, valuesy, 1)
    folded = np.append(folded, valuesz, 1)

    unfolded = initpoints

    return unfolded, folded
