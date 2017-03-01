import numpy as np
import copy

class Atom():
    element = ""
    fract = np.ndarray((1,3))
    cart = np.ndarray((1,3))

    def dirtocar(self,mat_dir_to_car):
        self.cart = np.array(np.inner(mat_dir_to_car,self.fract))

    def cartodir(self,mat_car_to_dir):
        self.fract = np.array(np.inner(mat_car_to_dir,self.cart))

    def move(self,vector):
        self.cart = self.cart + vector

    def remapID(self,oldID,newID):
        if (self.element == oldID):
            self.element = newID

class AtTraj():
    box = np.ndarray((3,3))
    abc = np.ndarray([1,3])
    ang = np.ndarray([1,3])
    description = ""
    atomlist = []
    mat_dir_to_car = np.zeros([3,3])
    mat_car_to_dir = np.zeros([3,3])

    def __init__(self):
        print 'Atomic trajectory initialized'


    def dirtocar(self):
        for atom in self.atomlist:
            atom.dirtocar(self.mat_dir_to_car)

    def cartodir(self):
        for atom in self.atomlist:
            atom.cartodir(self.mat_car_to_dir)

    def remapID(self,oldID,newID):
        for atom in self.atomlist:
            atom.remapID(oldID,newID)

    def vectortocar(self,inputvec):
        return np.inner(self.mat_dir_to_car,inputvec)

    def vectortodir(self,inputvec):
        return np.inner(self.mat_car_to_dir,inputvec)

    def move(self,vector):
        for atom in self.atomlist:
            atom.move(vector)

    def inbox(self):
        for atom in self.atomlist:
            atom.fract = atom.fract - np.floor(atom.fract)
            atom.dirtocar(self.mat_dir_to_car)


    def remove_overlap(self,cutoff):
        # Remove one of a pair of atoms that are within <cutoff> distance of each other
        # <cutoff> is given in Angstroms

        # Loop through the atomlist to identify atom pairs within a cutoff distance of each other and add one atom of the pair to a blacklist
        # Construct a whitelist of atoms in self.atomlist but not in the blacklist
        # Set self.atomlist equal to the whitelist

        num_of_initial_atoms = len(self.atomlist)
        removed_atoms = []

        for atom1 in self.atomlist:
            for atom2 in self.atomlist:
                if (atom1 not in removed_atoms) and (atom2 not in removed_atoms) and (atom1 != atom2) and (self.pbc_distance(atom1,atom2) < cutoff):
                    removed_atoms.append(atom2)

        num_of_removed_atoms = len(removed_atoms)
        remaining_atoms = [atom for atom in self.atomlist if atom not in removed_atoms]
        num_of_remaining_atoms = len(remaining_atoms)

        self.atomlist = remaining_atoms
        print num_of_removed_atoms,'atoms removed. Atomlist size reduced from',num_of_initial_atoms,'to',num_of_remaining_atoms



    def pbc_distance(self,atom1,atom2):
        atom1.fract = atom1.fract - np.floor(atom2.fract)
        atom2.fract = atom2.fract - np.floor(atom2.fract)
        diff0 = atom2.fract - atom1.fract
        diff1 = diff0 + 1
        diff2 = diff0 - 1
        diff0 = np.append([diff0],[diff1],axis=0)
        diff0 = np.append(diff0,[diff2],axis=0)
        diff = []
        diff.append(sorted(diff0[:,0],key = abs)[0])
        diff.append(sorted(diff0[:,1],key = abs)[0])
        diff.append(sorted(diff0[:,2],key = abs)[0])
        mindist_fract_vec = np.array(diff)
        mindist_cart_vec = np.dot(np.matrix.transpose(self.box),mindist_fract_vec)
        mindist = np.linalg.norm(mindist_cart_vec)
        return mindist



    def make_periodic(self,num_of_images):
        # Duplicate self.atomlist -- Dont want to (re)duplicate all the new atoms
        original_atomlist = list(self.atomlist)
        # Clear the original list
        self.atomlist = []
        # Loop over z,y and x respectively
        for zimage in range(0,num_of_images[2]):
            for yimage in range(0,num_of_images[1]):
                for ximage in range(0,num_of_images[0]):
                    for singleatom in original_atomlist:
                        newatom = copy.copy(singleatom)
                        newatom.fract = newatom.fract + np.array([ximage,yimage,zimage])
                        newatom.fract = np.divide(newatom.fract,num_of_images)
                        self.atomlist.append(newatom)

        # Adjust the box sizes and recompile transformation matrices
        self.box[0,:] = self.box[0,:] * num_of_images[0]
        self.box[1,:] = self.box[1,:] * num_of_images[1]
        self.box[2,:] = self.box[2,:] * num_of_images[2]
        self.make_dircar_matrices()
        self.dirtocar()

        # Cleanup
        del original_atomlist[:]



    def read_snapshot_lammps(self,filename):
        lammps_snapfile = open(filename,"r")
        lammps_snapfile.readline() # Comment line for Timestep
        self.stepcount = int(lammps_snapfile.readline().strip())

        lammps_snapfile.readline() # Comment line for Atom Count
        atomcount = int(lammps_snapfile.readline().strip())

        lammps_snapfile.readline() # Comment line for box bounds
        xlo_bound,xhi_bound,xy = map(float,lammps_snapfile.readline().split())
        ylo_bound,yhi_bound,xz = map(float,lammps_snapfile.readline().split())
        zlo_bound,zhi_bound,yz = map(float,lammps_snapfile.readline().split())

        xlo = xlo_bound - min(0.0,xy,xz,xy+xz)
        xhi = xhi_bound - max(0.0,xy,xz,xy+xz)
        ylo = ylo_bound - min(0.0,yz)
        yhi = yhi_bound - max(0.0,yz)
        zlo = zlo_bound
        zhi = zhi_bound

        self.box[0,:] = [xhi-xlo, 0.0, 0.0]
        self.box[1,:] = [xy, yhi-ylo, 0.0]
        self.box[2,:] = [xz, yz, zhi-zlo]

        self.make_dircar_matrices() # Uniform representation of box dimensions

        lammps_snapfile.readline() # Comment line for Atom positions
        for thisatomcount in range(0,atomcount):
            basisline = lammps_snapfile.readline()
            myatom = Atom()
            myatom.element, myatom.cart = [basisline.split()[0],  np.array(map(float,basisline.split()[1:4]))]
            self.atomlist.append(myatom)
        self.cartodir()

        lammps_snapfile.close()


    def read_snapshot_vasp(self,filename):
        vasp_snapfile = open(filename,"r")

        self.description = vasp_snapfile.readline().strip()
        mymultiplier = float(vasp_snapfile.readline())

        self.box[0,:] = map(float,vasp_snapfile.readline().split())
        self.box[1,:] = map(float,vasp_snapfile.readline().split())
        self.box[2,:] = map(float,vasp_snapfile.readline().split())
        self.box = self.box * mymultiplier

        self.make_dircar_matrices() # Uniform representation of box dimensions from POSCAR file

        basisline = vasp_snapfile.readline()
        atarray = basisline.split()

        basisline = vasp_snapfile.readline()
        atoms_of_type = map(int,basisline.split())

        isindirectcoords = vasp_snapfile.readline().lower().strip()[0]=='d'  # Check if the coordinates are in Direct or Cartesian

        if isindirectcoords:
            for index, numbers in enumerate(atoms_of_type):
                for thistype in range (0,numbers):
                    basisline = vasp_snapfile.readline()
                    myatom = Atom()
                    myatom.fract = np.array(map(float,basisline.split()))
                    myatom.element = atarray[index].upper()
                    self.atomlist.append(myatom)
            self.dirtocar() # Populate the cartesian position values from the fractional coordinates for each atom
        else:
            for index, numbers in enumerate(atoms_of_type):
                for thistype in range (0,numbers):
                    basisline = vasp_snapfile.readline()
                    myatom = Atom()
                    myatom.cart = np.array(map(float,basisline.split()))
                    myatom.element = atarray[index].upper()
                    self.atomlist.append(myatom)
            self.cartodir() # Populate fractional coordinates from the cartesian position of each atom

        vasp_snapfile.close()


    def make_dircar_matrices(self):
        # Convert box vectors into cell lengths (a,b,c) and angles (alpha, beta, gamma)
        a = np.linalg.norm(self.box[0,:])
        b = np.linalg.norm(self.box[1,:])
        c = np.linalg.norm(self.box[2,:])
        alpha = np.arccos( np.inner(self.box[1,:],self.box[2,:])/(b*c) )
        beta = np.arccos( np.inner(self.box[2,:],self.box[0,:])/(c*a) )
        gamma = np.arccos( np.inner(self.box[0,:],self.box[1,:])/(a*b) )

        self.abc = np.array([a, b, c])
        self.ang = np.array([alpha, beta , gamma])
        self.boxvolume = np.inner(self.box[0,:], np.cross(self.box[1,:],self.box[2,:]))

        self.mat_dir_to_car[0,:] = [a , b * np.cos(gamma) , c * np.cos(beta)]
        self.mat_dir_to_car[1,:] = [0 , b * np.sin(gamma) , c * ((np.cos(alpha) - (np.cos(beta)*np.cos(gamma)))/ np.sin(gamma))   ]
        self.mat_dir_to_car[2,:] = [0 , 0 , self.boxvolume / (a * b * np.sin(gamma))   ]

        self.mat_car_to_dir[0,:] = [1.0/a ,  0.0 - (np.cos(gamma)/(a * np.sin(gamma))) , b * c * ((np.cos(alpha)*np.cos(gamma)) - np.cos(beta))/(self.boxvolume * np.sin(gamma))       ]
        self.mat_car_to_dir[1,:] = [0.0 ,  1.0 / (b * np.sin(gamma)) , a * c * ((np.cos(beta)*np.cos(gamma)) - np.cos(alpha))/(self.boxvolume * np.sin(gamma))       ]
        self.mat_car_to_dir[2,:] = [0.0 ,  0.0 , a * b * np.sin(gamma) / self.boxvolume]


    def sort_by_element(self):
        self.atomlist.sort(key = lambda x: x.element)


    def write_snapshot_vasp(self,filename,write_in_direct):
        vasp_snapfile = open(filename,"w")
        vasp_snapfile.write(self.description+'\n')
        vasp_snapfile.write('1.000000\n') # Default multiplier for all VASP files
        np.savetxt(vasp_snapfile,self.box,fmt='%19.16f', delimiter = "   ", newline = "\n")

        # Sort atoms by element before counting number of atoms by element
        self.atomlist.sort(key = lambda x: x.element)
        uniquesdict = {}
        for singleatom in self.atomlist:
            uniquesdict[singleatom.element] = singleatom.element
        uniqueslist = uniquesdict.values()

        uniquesdict = {}
        uniquesdict1 = {}
        for uniqueelement in uniqueslist:
            uniquesdict[uniqueelement] = str(len([p for p in self.atomlist if p.element == uniqueelement]))
            uniquesdict1[uniqueelement] = uniqueelement
        uniqueslist = uniquesdict.values()
        uniqueslist1 = uniquesdict1.values()

        vasp_snapfile.write("  ".join(uniqueslist1).title()+'\n')
        vasp_snapfile.write("  ".join(uniqueslist)+'\n')

        if write_in_direct:
            vasp_snapfile.write('Direct\n')
            for uniqueelement in uniqueslist1:
                subsetofatomlist = (atoms for atoms in self.atomlist if atoms.element == uniqueelement)
                for singleatom in subsetofatomlist:
                    np.savetxt(vasp_snapfile,singleatom.fract[None],fmt='%19.16f', delimiter = "   ", newline = "\n ") # np.savetxt has problems with 1D array writing
        else:
            vasp_snapfile.write('Cartesian\n')
            for uniqueelement in uniqueslist1:
                subsetofatomlist = (atoms for atoms in self.atomlist if atoms.element == uniqueelement)
                for singleatom in subsetofatomlist:
                    np.savetxt(vasp_snapfile,singleatom.cart[None],fmt='%19.16f', delimiter = "   ", newline = "\n ") # np.savetxt has problems with 1D array writing



        vasp_snapfile.close()





def isSierpinskiCarpetFilled(level,coords):
    # Calculate if the given fractional coordinates correspond to a filed pixel in the Sierpinksi carpet fractal of a given level

    # Multiply the fractional coordinate with 3^n (for n = level..1) and check to see if it leaves a reminder of 1 upon division by 3 (i.e. it is the middle cell at any level)

    multiplier = 3**level

    x = int(coords[0]*multiplier)
    y = int(coords[1]*multiplier)

    while True:
        if (x==0 and y==0):
            break

        if (x%3==1 and y%3==1):
            return False

        x = int(x/3)
        y = int(y/3)

    return True
