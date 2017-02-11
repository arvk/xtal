import numpy as np

class AtTraj:
    # Define variables
    box = np.ndarray((3,3))
    basisa = box[0,:]
    basisb = box[1,:]
    basisc = box[2,:]

    def __init__(self):
        print 'Atomic trajectory initialized'

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
                myatom.dirtocar(self.box)
                self.atomlist.append(myatom)

        vasp_snapfile.close()


    def sort_by_element(self):
        self.atomlist.sort(key = lambda x: x.element)

