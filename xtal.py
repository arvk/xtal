import numpy as np

class Atom():
    afrac, bfrac, cfrac = [0.0, 0.0, 0.0]
    xpos, ypos, zpos = [0.0, 0.0, 0.0]
    element = ""

    def __init__(self):
        print 'New atom'

    def dirtocar(self,box):
        self.xpos = self.afrac*box[0,0] + self.bfrac*box[1,0] + self.cfrac*box[2,0]
        self.ypos = self.afrac*box[0,1] + self.bfrac*box[1,1] + self.cfrac*box[2,1]
        self.zpos = self.afrac*box[0,2] + self.bfrac*box[1,2] + self.cfrac*box[2,2]



class AtTraj():
    box = np.ndarray((3,3))
    description = ""
    atomlist = []

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

        self.boxvolume = np.inner(self.box[0,:], np.cross(self.box[1,:],self.box[2,:]))


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
