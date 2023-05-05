from periodictable import core, mass, density
import copy
import numpy as np
from tess import Container
import xtal

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


def make_polycrystal(grain_snapshot, num_grains, box_size, remove_duplicates = None):
    '''Make polycrystalline structure out of a given grain atomic structure'''
    # Identify voronoi points uniformly distributed within the desired box size
    random_positions = np.random.rand(num_grains,3)
    initial_voronoi_points = [position*box_size.tolist() for position in random_positions]

    # Generate the output trajectory
    v = xtal.AtTraj()
    v.create_snapshot(xtal.Snapshot)
    v.abc = box_size
    v.ang = (np.pi/2.0) * np.array([1.0,1.0,1.0])
    v.abc_to_box()
    v.make_dircar_matrices()

    # Generate Voronoi structure using tess
    cntr = Container(initial_voronoi_points, limits=tuple(box_size), periodic=True)
    cell_volumes = [cell.volume() for cell in cntr]

    # Replicate grain snapshot to fill the largest voronoi volume
    max_vol = np.max(cell_volumes)
    unit_cell_volume = grain_snapshot.trajectory.boxvolume
    volume_ratio = max_vol / unit_cell_volume
    rep_ratio = np.cbrt(volume_ratio)
    rep_ratio = int(np.ceil(rep_ratio)) * 3
    grain_snapshot.trajectory.make_periodic([rep_ratio, rep_ratio, rep_ratio])

    grain_snapshot.dirtocar()
    #print(grain_snapshot.atomlist[0].cart)
    #print(grain_snapshot.atomlist[1].cart)

    # Identify center of unit cell
    uc_center = np.array([0.0,0.0,0.0])
    for atom in grain_snapshot.atomlist:
        uc_center += atom.cart
    uc_center /= len(grain_snapshot.atomlist)

    # Rotate unit cell and template inside each voronoi volume
    for cellID, cell in enumerate(cntr):
        grain_snapshot.dirtocar()
        angle_z = np.random.random() * np.pi
        angle_y = np.random.random() * np.pi
        angle_x = np.random.random() * np.pi
        grain_snapshot.rotate_euler(uc_center, angle_z, angle_y, angle_x)
        #centroid = cntr[cellID].pos
        centroid = cntr[cellID].centroid()
        grain_snapshot.move(centroid - uc_center)

        for atom in grain_snapshot.atomlist:
            same_side = True
            for face in cntr[cellID].face_vertices():
                pt1 = cntr[cellID].vertices()[face[0]]
                pt2 = cntr[cellID].vertices()[face[1]]
                pt3 = cntr[cellID].vertices()[face[2]]
                offset_pt2 = np.array(pt2) - np.array(pt1)
                offset_pt3 = np.array(pt3) - np.array(pt1)
                plane = np.cross(offset_pt2, offset_pt3)
                side_centroid = np.sign(np.dot(plane, centroid) - np.dot(plane, pt1))
                side_atom = np.sign(np.dot(plane, atom.cart) - np.dot(plane, pt1))
                if side_centroid != side_atom:
                    same_side = False # This atom is not inside the Voronoi volume
                    break

            if same_side:
                vatom = v.snaplist[0].create_atom(xtal.Atom)
                vatom.element = atom.element
                vatom.cart = atom.cart
                vatom.cartodir()

    # Wrap atoms within box and remove near-overlapping atoms near the grain boundaries
    v.inbox()
    if remove_duplicates:  # If a remove_duplicates distance is provided
        v.snaplist[0].remove_duplicates(remove_duplicates)
    return v
