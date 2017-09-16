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

    def test_remap_id(self):
        """Test if VASP5 POSCARs can be read"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.remap_id('S','TE')
        assert len([atom for atom in u.snaplist[0].atomlist if atom.element == 'S']) == 0
        assert len([atom for atom in u.snaplist[0].atomlist if atom.element == 'TE']) == 2
