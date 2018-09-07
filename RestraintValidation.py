#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import ntpath
import numpy

PY3 = (sys.version_info[0] == 3)

(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))

try:
    import pynmrsta

except ImportError as e:
    sys.path.append(scriptPath + '/PyNMRSTAR')
    try:
        import bmrb as pynmrstar

        print("Using local pynmrstar")
    except ImportError as e:
        print(scriptPath, scriptName)
        print("ERROR: PyNMRSTAR parser is not available")
        print(str(e))
        exit(1)

try:
    from mmcif.io.PdbxReader import PdbxReader
except ImportError as e:
    sys.path.append(scriptPath + '/py-mmcif')
    try:
        from mmcif.io.PdbxReader import PdbxReader

        print("Using local mmcif parser")
    except ImportError as e:
        print(scriptPath, scriptName)
        print("ERROR: mmcif parser is not available")
        print(str(e))
        exit(1)

__version__ = "v1.0"


class RestraintValidation:

    def __init__(self):
        print (pynmrstar.__version__)
        self.cfile = 'test/data/1nk2.cif'
        self.rfile = 'test/data/1nk2.str'
        self.cif_data = None
        self.max_models = None
        self.read_ciffile()
        self.star_data = None
        self.rest_flag = False
        self.read_starfile()
        self.ang_rest = None
        self.get_angle_restraints()



    @staticmethod
    def r6sum(dist_list):
        return (sum([i ** (-6.) for i in dist_list])) ** (-1. / 6.)

    @staticmethod
    def r6average(dist_list):
        return (sum([i ** (-6.) for i in dist_list]) / float(len(dist_list))) ** (-1. / 6.)

    def check_files(self):
        """ Check the existence of input files"""
        if os.path.isfile(self.cfile):
            print('Coordinate file found {}'.format(self.cfile))
        else:
            print('Cooriantefile NOT found')
            return False
        if os.path.isfile(self.rfile):
            print('Restraint file found {}'.format(self.rfile))
        else:
            print('Restraint file NOT found')
            return False
        return True

    def read_ciffile(self):
        """ Reads the input coordinate CIF file using pdbx lightweight parser
        ref: http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html """
        self.cif_data = []
        ifh = open(self.cfile, 'r')
        pRd = PdbxReader(ifh)
        pRd.read(self.cif_data)
        ifh.close()
        c0 = self.cif_data[0]
        atom_site = c0.getObj('atom_site')
        self.max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
        if self.max_models == 1:
            print('Coordinate file has only one model')
        elif self.max_models == 0:
            print('Coordinate file has zero models')
        else:
            print('Coordinate file has {} models'.format(self.max_models))

    def read_starfile(self):
        """ Reads the input NMR-STAR file using pynmrstar library ref : https://github.com/uwbmrb/PyNMRSTAR """
        self.star_data = pynmrstar.Entry.from_file(self.rfile)
        cat_list = [saveframe.category for saveframe in self.star_data]
        rest_info = [False, False]
        if 'general_distance_constraints' in cat_list:
            print ('general_distance_constraints saveframe found')
            rest_info[0] = True
        else:
            print ('general_distance_constraints saveframe not found')
        if 'torsion_angle_constraints' in cat_list:
            print ('torsion_angle_constraints saveframe found')
            rest_info[1] = True
        else:
            print ('torsion_angle_constraints saveframe not found')

        if 'general_distance_constraints' not in cat_list and 'torsion_angle_constraints' not in cat_list:
            print (
                'Both general_distance_constraints and torsion_angle_constraints are missing in the STAR file')
            rest_info = [False, False]
        self.rest_flag = rest_info

    def get_seq_len(self):
        seq = self.star_data.get_tag('_Chem_comp_assembly.Comp_index_ID')
        return len(seq)

    def get_coordinates(self, modelID):

        """Creates a dictionary of coordinates gor a given modelID.
        co[atom] = array([x,y,z])
        atoms are identified by (label_seq_id,label_entity_id,label_comp_id,lable_atom_id)"""

        c0 = self.cif_data[0]
        atom_site = c0.getObj('atom_site')
        colnames = atom_site.getAttributeList()
        modelid = colnames.index('pdbx_PDB_model_num')
        self.max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
        xid = colnames.index('Cartn_x')
        yid = colnames.index('Cartn_y')
        zid = colnames.index('Cartn_z')
        atomid = colnames.index('label_atom_id')
        compid = colnames.index('label_comp_id')
        asymid = colnames.index('label_asym_id')
        entityid = colnames.index(
            'label_entity_id')  # asymid is used instead of entity id because of different definition of cif and star
        seqid = colnames.index('label_seq_id')
        co = {}
        for dat in atom_site.getRowList():
            if int(dat[modelid]) == modelID:
                co[(dat[seqid], dat[asymid], dat[compid], dat[atomid])] = numpy.array(
                    [float(dat[xid]), float(dat[yid]), float(dat[zid])])
        return co

    @staticmethod
    def get_dihedral_angle(c1, c2, c3, c4):

        """ Calculates the dihedral angle from the given four coordinate values.
        Each coordinate is an array of x,y,z. Returns angle in degrees"""

        bv12 = c1 - c2
        bv32 = c3 - c2
        bv43 = c4 - c3
        pv13 = numpy.cross(bv12, bv32)
        pv24 = numpy.cross(bv43, bv32)
        pro = numpy.dot(pv13, pv24)
        sqdist13 = numpy.dot(pv13, pv13)
        sqdist24 = numpy.dot(pv24, pv24)
        cosin = pro / sqrt(sqdist13 * sqdist24)
        cosin - min(1.0, max(-1.0, cosin))
        angle = acos(cosin)

        if numpy.dot(pv13, numpy.cross(pv24, bv32)) < 0:
            angle = -angle
        return round(numpy.degrees(angle), 4)

    @staticmethod
    def get_distance(c1, c2):

        """ Calculates the distance between two coordinate values.
        Each coordinate is an array of x,y,z. Returns distance in A assuming the input coordinates are in A"""

        return numpy.linalg.norm(c1 - c2)

    def get_angle_restraints(self):

        """Creates a dictionary of angle restraints from _Torsion_angle_constraint loop.
        ang_rest[atom 1,atom 2, atom 3, atom 4] = [lower bound, upper bound]
        atoms are identified by (Comp_index_ID,Entity_assembly_ID,Comp_ID,Atom_ID)"""

        ang = self.star_data.get_loops_by_category('_Torsion_angle_constraint')
        self.ang_rest = []
        for dl in ang:
            self.ang_rest.append({})
            restid = dl.get_tag_names().index("_Torsion_angle_constraint.ID")
            rest_name = dl.get_tag_names().index("_Torsion_angle_constraint.Torsion_angle_name")
            atomid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_1")
            atomid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_2")
            atomid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_3")
            atomid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Atom_ID_4")
            compid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_1")
            compid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_2")
            compid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_3")
            compid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_ID_4")
            # entityid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_1")
            # entityid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_2")
            # entityid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_3")
            # entityid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Entity_assembly_ID_4")
            entityid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_1")
            entityid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_2")
            entityid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_3")
            entityid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Auth_asym_ID_4")
            seqid1 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_1")
            seqid2 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_2")
            seqid3 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_3")
            seqid4 = dl.get_tag_names().index("_Torsion_angle_constraint.Comp_index_ID_4")
            lb = dl.get_tag_names().index("_Torsion_angle_constraint.Angle_lower_bound_val")
            ub = dl.get_tag_names().index("_Torsion_angle_constraint.Angle_upper_bound_val")
            for i in dl:
                try:
                    lbv = float(i[lb])
                except ValueError:
                    lbv = -999
                try:
                    ubv = float(i[ub])
                except ValueError:
                    ubv = -999
                self.ang_rest[-1][int(i[restid])] = [i[rest_name], (i[seqid1], i[entityid1], i[compid1], i[atomid1]),
                                                     (i[seqid2], i[entityid2], i[compid2], i[atomid2]),
                                                     (i[seqid3], i[entityid3], i[compid3], i[atomid3]),
                                                     (i[seqid4], i[entityid4], i[compid4], i[atomid4]), (lbv, ubv)]
        print ('Number of angle restraints : {}'.format([len(i) for i in self.ang_rest]))


if __name__ == "__main__":
    p = RestraintValidation()
