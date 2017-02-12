import numpy as np

class Atom():
    element = ""
    fract = np.ndarray((1,3))
    cart = np.ndarray((1,3))

    def __init__(self):
        print 'New atom'

    def dirtocar(self,mat_dir_to_car):
        self.cart = np.array(np.inner(mat_dir_to_car,self.fract))

    def cartodir(self,mat_car_to_dir):
        self.fract = np.array(np.inner(mat_car_to_dir,self.cart))

    def move(self,vector):
        self.cart = self.cart + vector

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

    def vectortocar(self,inputvec):
        return np.inner(self.mat_dir_to_car,inputvec)

    def vectortodir(self,inputvec):
        return np.inner(self.mat_car_to_dir,inputvec)

    def move(self,vector):
        for atom in self.atomlist:
            atom.move(vector)


    def read_snapshot_vasp(self,filename):
        vasp_snapfile = open(filename,"r")

        self.description = vasp_snapfile.readline().strip()
        mymultiplier = float(vasp_snapfile.readline())

        self.box[0,:] = map(float,vasp_snapfile.readline().split())
        self.box[1,:] = map(float,vasp_snapfile.readline().split())
        self.box[2,:] = map(float,vasp_snapfile.readline().split())
        self.box = self.box * mymultiplier

        self.boxvolume = np.inner(self.box[0,:], np.cross(self.box[1,:],self.box[2,:]))

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
                subsetofatomlist = (atoms for atoms in self.atomlist if atoms.element == uniqueelement.upper())
                for singleatom in subsetofatomlist:
                    np.savetxt(vasp_snapfile,singleatom.cart[None],fmt='%19.16f', delimiter = "   ", newline = "\n ") # np.savetxt has problems with 1D array writing



        vasp_snapfile.close()
