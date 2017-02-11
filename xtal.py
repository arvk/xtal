import numpy as np

class Atom():
    afrac, bfrac, cfrac = [0.0, 0.0, 0.0]
    xpos, ypos, zpos = [0.0, 0.0, 0.0]
    element = ""

    def __init__(self):
        print 'New atom'

    def dirtocar(self,mat_dir_to_car):
        [self.xpos, self.ypos, self.zpos] = np.inner(mat_dir_to_car,[self.afrac, self.bfrac, self.cfrac])

    def cartodir(self,mat_car_to_dir):
        [self.afrac, self.bfrac, self.cfrac] = np.inner(mat_car_to_dir,[self.xpos, self.ypos, self.zpos])


class AtTraj():
    box = np.ndarray((3,3))
    description = ""
    atomlist = []

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



    def read_snapshot_vasp(self,filename):
        vasp_snapfile = open(filename,"r")

        self.description = vasp_snapfile.readline().strip()
        mymultiplier = float(vasp_snapfile.readline())

        basisline = vasp_snapfile.readline()
        self.box[0,0], self.box[0,1], self.box[0,2] = map(float,basisline.split())
        basisline = vasp_snapfile.readline()
        self.box[1,0], self.box[1,1], self.box[1,2] = map(float,basisline.split())
        basisline = vasp_snapfile.readline()
        self.box[2,0], self.box[2,1], self.box[2,2] = map(float,basisline.split())
        self.box = self.box * mymultiplier

        self.boxvolume = np.inner(self.box[0,:], np.cross(self.box[1,:],self.box[2,:]))

        self.make_dircar_matrices() # Uniform representation of box dimensions from POSCAR file

        basisline = vasp_snapfile.readline()
        atarray = basisline.split()

        basisline = vasp_snapfile.readline()
        atoms_of_type = map(int,basisline.split())

        isindirectcoords = vasp_snapfile.readline()
        isindirectcoords = isindirectcoords[0].lower() == 'd'

        for index, numbers in enumerate(atoms_of_type):
            for thistype in range (0,numbers):
                basisline = vasp_snapfile.readline()
                myatom = Atom()
                myatom.afrac, myatom.bfrac, myatom.cfrac = map(float,basisline.split())
                myatom.element = atarray[index].upper()
                self.atomlist.append(myatom)

        self.dirtocar() # Populate the cartesian position values from the fractional coordinates for each atom

        vasp_snapfile.close()


    def make_dircar_matrices(self):
        # Convert cell vectors a, b, c into cell lengths and angles
        self.boxa = np.linalg.norm(self.box[0,:])
        self.boxb = np.linalg.norm(self.box[1,:])
        self.boxc = np.linalg.norm(self.box[2,:])
        self.boxalpha = np.arccos( np.inner(self.box[1,:],self.box[2,:])/(self.boxb*self.boxc) )
        self.boxbeta = np.arccos( np.inner(self.box[2,:],self.box[0,:])/(self.boxc*self.boxa) )
        self.boxgamma = np.arccos( np.inner(self.box[0,:],self.box[1,:])/(self.boxa*self.boxb) )

        self.mat_dir_to_car = np.zeros([3,3])
        self.mat_dir_to_car[0,:] = [self.boxa , self.boxb * np.cos(self.boxgamma) , self.boxc * np.cos(self.boxbeta)]
        self.mat_dir_to_car[1,:] = [0 , self.boxb * np.sin(self.boxgamma) , self.boxc * ((np.cos(self.boxalpha) - (np.cos(self.boxbeta)*np.cos(self.boxgamma)))/ np.sin(self.boxgamma))   ]
        self.mat_dir_to_car[2,:] = [0 , 0 , self.boxvolume / (self.boxa * self.boxb * np.sin(self.boxgamma))   ]

        self.mat_car_to_dir = np.zeros([3,3])
        self.mat_car_to_dir[0,:] = [1.0/self.boxa ,  0.0 - (np.cos(self.boxgamma)/(self.boxa * np.sin(self.boxgamma))) , self.boxb * self.boxc * ((np.cos(self.boxalpha)*np.cos(self.boxgamma)) - np.cos(self.boxbeta))/(self.boxvolume * np.sin(self.boxgamma))       ]
        self.mat_car_to_dir[1,:] = [0.0 ,  1.0 / (self.boxb * np.sin(self.boxgamma)) , self.boxa * self.boxc * ((np.cos(self.boxbeta)*np.cos(self.boxgamma)) - np.cos(self.boxalpha))/(self.boxvolume * np.sin(self.boxgamma))       ]
        self.mat_car_to_dir[2,:] = [0.0 ,  0.0 , self.boxa * self.boxb * np.sin(self.boxgamma) / self.boxvolume]


    def sort_by_element(self):
        self.atomlist.sort(key = lambda x: x.element)


    def write_snapshot_vasp(self,filename,write_in_direct):
        vasp_snapfile = open(filename,"w")
        vasp_snapfile.write(self.description+'\n')
        vasp_snapfile.write('1.000000\n') # Default multiplier for all VASP files
        np.savetxt(vasp_snapfile,self.box,fmt='%19.16f', delimiter = "   ", newline = "\n")

        # Sort atoms by element before counting number of atoms by element
        #self.atomlist.sort_by_element()
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
                    vasp_snapfile.write('{0}  {1}  {2}\n'.format(singleatom.afrac,singleatom.bfrac,singleatom.cfrac))
        else:
            vasp_snapfile.write('Cartesian\n')
            for uniqueelement in uniqueslist1:
                subsetofatomlist = (atoms for atoms in self.atomlist if atoms.element == uniqueelement.upper())
                for singleatom in subsetofatomlist:
                    vasp_snapfile.write('{0}  {1}  {2}\n'.format(singleatom.xpos,singleatom.ypos,singleatom.zpos))



        vasp_snapfile.close()
