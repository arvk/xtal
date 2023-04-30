import xtal
import numpy as np
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
        v = xtal.AtTraj()
        v.read_snapshot_vasp('tests/POSCAR')
        u_mo_atoms = [atom for atom in u.snaplist[0].atomlist if atom.element == "Mo"]
        v_mo_atoms = [atom for atom in v.snaplist[0].atomlist if atom.element == "Mo"]
        assert np.linalg.norm(u_mo_atoms[0].cart - v_mo_atoms[0].cart) < 0.0001
        os.remove('tests/POSCAR')

    def test_write_vasp_poscar_direct(self):
        """Test if VASP5 POSCARs can be written in fractional coordinates"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.cartesian.unitcell')
        u.snaplist[0].write_snapshot_vasp('tests/POSCAR',True)
        v = xtal.AtTraj()
        v.read_snapshot_vasp('tests/POSCAR')
        u_mo_atoms = [atom for atom in u.snaplist[0].atomlist if atom.element == "Mo"]
        v_mo_atoms = [atom for atom in v.snaplist[0].atomlist if atom.element == "Mo"]
        assert np.linalg.norm(u_mo_atoms[0].fract - v_mo_atoms[0].fract) < 0.0001
        os.remove('tests/POSCAR')


# Testing PWP integration
class TestPWP(object):

    def test_read_pwp_trajectory(self):
        """Test if a PWP trajectory can be read"""
        u = xtal.AtTraj()
        u.read_trajectory_pwp('tests/pwp_aimd')
        num_snapshots = len(u.snaplist)
        num_atoms_last_snapshot = len(u.snaplist[-1].atomlist)
        assert (num_snapshots, num_atoms_last_snapshot) == (101, 54)


# Testing General trajectory methods
class TestGeneral(object):

    def test_remap_id(self):
        """Test if trajectory elements can be remapped"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.remap_id('S','Te')
        no_S_atoms = len([atom for atom in u.snaplist[0].atomlist if atom.element == 'S'])
        no_Te_atoms = len([atom for atom in u.snaplist[0].atomlist if atom.element == 'Te'])
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

    def test_remove_duplicates(self):
        """Test if duplicate atoms are removed based on provided cutoff"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.duplicates.unitcell')
        snapshot = u.snaplist[0]
        raw_num_atoms = len(snapshot.atomlist)

        # Remove all duplicate atoms closer than 0.11
        snapshot.remove_duplicates(0.11)
        pt1_pass1_num_atoms = len(snapshot.atomlist)
        snapshot.remove_duplicates(0.11) # second pass shold not change anything
        pt1_pass2_num_atoms = len(snapshot.atomlist)

        # Remove all duplicate atoms closer than 0.21
        snapshot.remove_duplicates(0.21)
        pt2_pass1_num_atoms = len(snapshot.atomlist)
        snapshot.remove_duplicates(0.21) # second pass shold not change anything
        pt2_pass2_num_atoms = len(snapshot.atomlist)

        # Remove all duplicate atoms closer than 0.31
        snapshot.remove_duplicates(0.31)
        pt3_pass1_num_atoms = len(snapshot.atomlist)
        snapshot.remove_duplicates(0.31) # second pass shold not change anything
        pt3_pass2_num_atoms = len(snapshot.atomlist)

        assert (raw_num_atoms, pt1_pass1_num_atoms, pt1_pass2_num_atoms, \
                pt2_pass1_num_atoms, pt2_pass2_num_atoms, \
                pt3_pass1_num_atoms, pt3_pass2_num_atoms) == (6, 5, 5, 4, 4, 3, 3)

    def test_move_atoms(self):
        """Test if overlapping atoms are removed based on provided cutoff"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.dirtocar()
        snapshot = u.snaplist[0]
        mo_atom = [atom for atom in snapshot.atomlist if atom.element == "Mo"][0]
        s_atom_top = [atom for atom in snapshot.atomlist if atom.fract[2] > 0.15][0]
        s_atom_bottom = [atom for atom in snapshot.atomlist if atom.fract[2] < 0.1][0]

        # Establish baselines
        top_dist = snapshot.pbc_distance(mo_atom,s_atom_top)
        bot_dist = snapshot.pbc_distance(mo_atom,s_atom_bottom)

        # Move the top and bottom S atoms by different distances. This should increase pbc_distances
        s_atom_top.move(np.array([0.0,0.0,5.0]))
        s_atom_bottom.move(np.array([0.0,0.0,-6.0]))
        u.cartodir()
        m1_top_dist = snapshot.pbc_distance(mo_atom,s_atom_top)
        m1_bot_dist = snapshot.pbc_distance(mo_atom,s_atom_bottom)

        # Move the entire trajectory. This should not change pbc_distances
        u.move(np.array([0.0,0.0,9.0]))
        u.cartodir()
        m2_top_dist = snapshot.pbc_distance(mo_atom,s_atom_top)
        m2_bot_dist = snapshot.pbc_distance(mo_atom,s_atom_bottom)

        # Test all distances. Initially 2.4, Expands beyond 5.0. And then stays the same
        assert (top_dist < 2.41 and bot_dist < 2.41 and \
                m1_top_dist > 5.0 and m1_bot_dist > 6.0 and \
                np.isclose(m2_top_dist, m1_top_dist, atol=1e-5) and \
                np.isclose(m2_bot_dist, m1_bot_dist, atol=1e-5))


    def test_inbox(self):
        """Test if atoms can be folded back into the supercell"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.cartesian.unitcell')
        u.snaplist[0].atomlist[0].move(np.array([0.0,0.0,30.0]))
        u.cartodir()
        u.inbox()
        zpos_fract = u.snaplist[0].atomlist[0].fract[2]
        zpos_cart = u.snaplist[0].atomlist[0].cart[2]
        assert (zpos_fract < 0.3 and zpos_cart < 10.0)


    def test_boxtoabc(self):
        """Test if atoms can be folded back into the supercell"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.box_to_abc()
        assert np.allclose(u.abc, np.array([3.16, 3.16, 27.16]), atol=0.1) and \
               np.allclose(u.ang, np.array([np.pi/2, np.pi/2, 2*np.pi/3]), atol=0.1)

    def test_fractal(self):
        """Test fractal method"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.make_periodic(np.array([100,100,1]))
        lv0 = len(u.snaplist[0].atomlist)
        selatoms = [atom for atom in u.snaplist[0].atomlist if xtal.is_sierpinski_carpet_filled(1,atom.fract)]
        lv1 = len(selatoms)
        selatoms = [atom for atom in u.snaplist[0].atomlist if xtal.is_sierpinski_carpet_filled(2,atom.fract)]
        lv2 = len(selatoms)
        selatoms = [atom for atom in u.snaplist[0].atomlist if xtal.is_sierpinski_carpet_filled(3,atom.fract)]
        lv3 = len(selatoms)
        selatoms = [atom for atom in u.snaplist[0].atomlist if xtal.is_sierpinski_carpet_filled(4,atom.fract)]
        lv4 = len(selatoms)

        assert (lv0, lv1, lv2, lv3, lv4) == (30000, 26634, 23730, 21210, 18813)


    def test_identisort_cartesian(self):
        """Test non periodic identisort"""
        u = xtal.AtTraj()
        u.read_snapshot_vasp('tests/POSCAR.VASP5.unitcell')
        u.read_snapshot_vasp('tests/POSCAR.VASP5.randomized.unitcell')
        u.dirtocar()

        snap0 = u.snaplist[0]
        snap1 = u.snaplist[1]

        u.identical_sort(snap0, snap1, periodic=False)
        assert np.linalg.norm(snap0.atomlist[0].cart - snap1.atomlist[0].cart) < 0.01 and \
               np.linalg.norm(snap0.atomlist[1].cart - snap1.atomlist[1].cart) < 0.01 and \
               np.linalg.norm(snap0.atomlist[2].cart - snap1.atomlist[2].cart) < 0.01
