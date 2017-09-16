import xtal
import numpy as np

# Testing VASP integration
class TestVASP(object):

    def test_read_vasp_poscar(self):
        """Test if VASP5 POSCARs can be read"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        assert len(u.snaplist[0].atomlist) == 3

    def test_make_periodic_vasp_poscar(self):
        """Test if VASP5 POSCARs an be PBC replicated"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.make_periodic(np.array([5, 5, 3]))
        assert len(u.snaplist[0].atomlist) == 225



# Testing General trajectory methods
class TestGeneral(object):

    def test_remap_id(self):
        """Test if VASP5 POSCARs can be read"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.remap_id('S','TE')
        no_S_atoms = len([atom for atom in u.snaplist[0].atomlist if atom.element == 'S'])
        no_Te_atoms = len([atom for atom in u.snaplist[0].atomlist if atom.element == 'TE'])
        assert (no_S_atoms, no_Te_atoms) == (0,2)


    def test_pbc_distance(self):
        """Test if distances between atoms are calculated correctly"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        mo_atom = u.snaplist[0].atomlist[0]
        s_atom = u.snaplist[0].atomlist[1]
        same_atom_distance = u.snaplist[0].pbc_distance(mo_atom,mo_atom)
        different_atom_distance = u.snaplist[0].pbc_distance(mo_atom,s_atom)
        assert (same_atom_distance, different_atom_distance) == (0.0, 2.4084020256748055)
