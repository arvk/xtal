import xtal
import numpy as np
import filecmp
import os

# Testing VASP integration
class TestVASP(object):

    def test_read_vasp_poscar_direct(self):
        """Test if VASP5 POSCARs in direct coordinates can be read"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        assert len(u.snaplist[0].atomlist) == 3

    def test_read_vasp_poscar_cartesian(self):
        """Test if VASP5 POSCARs in cartesian coordinates can be read"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.cartesian.unitcell')
        assert len(u.snaplist[0].atomlist) == 3

    def test_make_periodic_vasp_poscar(self):
        """Test if VASP5 POSCARs an be PBC replicated"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.make_periodic(np.array([5, 5, 3]))
        assert len(u.snaplist[0].atomlist) == 225

    def test_write_vasp_poscar_cartesian(self):
        """Test if VASP5 POSCARs can be written in cartesian coordinates"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.snaplist[0].write_snapshot_vasp('tests/POSCAR',False)
        assert filecmp.cmp('tests/POSCAR','tests/POSCAR.VASP5.cartesian.unitcell')
        os.remove('tests/POSCAR')


# Testing General trajectory methods
class TestGeneral(object):

    def test_remap_id(self):
        """Test if trajectory elements can be remapped"""
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
        assert (same_atom_distance == 0.0 and different_atom_distance > 2.40 and different_atom_distance < 2.41)

    def test_dirtocar(self):
        """Test if fractional units can be converted to cartesian coordinates"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.dirtocar()
        zpos_atom_one = u.snaplist[0].atomlist[0].cart[2]
        assert (zpos_atom_one > 3.0 and zpos_atom_one < 3.1)
