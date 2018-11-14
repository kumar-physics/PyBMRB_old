#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import ntpath
import numpy
from math import sqrt,acos,pi
from cmath import exp
import operator

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


class Cyrange:

    def __init__(self,ciffile):
        self.cfile = ciffile
        self.read_ciffile()
        self.dihedral_angles = {
            'ALA':[['C','N','CA','C'],['N','CA','C','N']],
            'ARG':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','NE'],['CG','CD','NE','CZ'],['CD','NE','CZ','NH1']],
            'ASN':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','OD1']],
            'ASP':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','OD1']],
            'CYS':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','SG']],
            'GLN':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','OE1']],
            'GLU':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','OE1']],
            'GLY':[['C','N','CA','C'],['N','CA','C','N']],
            'HIS':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','ND1']],
            'ILE':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG1'],['CA','CB','CG1','CD1']],
            'LEU':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD1']],
            'LYS':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD'],['CB','CG','CD','CE'],['CG','CD','CE','NZ']],
            'MET':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','SD'],['CB','CG','SD','CE']],
            'PHE':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD1']],
            'PRO':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD']],
            'SER':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','OG']],
            'THR':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','OG1']],
            'TRP':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD1']],
            'TYR':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG'],['CA','CB','CG','CD1']],
            'VAL':[['C','N','CA','C'],['N','CA','C','N'],['N','CA','CB','CG1']],
        }
        self.angle_order_parameter()



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

    def get_seq_from_cif(self):
        self.cif_data = []
        ifh = open(self.cfile, 'r')
        pRd = PdbxReader(ifh)
        pRd.read(self.cif_data)
        ifh.close()
        c0 = self.cif_data[0]
        seq_dat = c0.getObj('entity_poly_seq')
        colnames = seq_dat.getAttributeList()
        entity_id = colnames.index('entity_id')
        seq_id = colnames.index('num')
        res_id = colnames.index('mon_id')
        max_entity = int(seq_dat.getValue('entity_id',-1))
        seq_dict=[]
        for i in range(max_entity):
            seq_dict.append({})
        for dat in seq_dat.getRowList():
            seq_dict[int(dat[entity_id])-1][int(dat[seq_id])]=dat[res_id]
        ent_dat = c0.getObj('entity_poly')
        colnames = ent_dat.getAttributeList()
        entity_id = colnames.index('entity_id')
        typ = colnames.index('type')
        entity_dict = {}
        for dat in ent_dat.getRowList():
            entity_dict[int(dat[entity_id])] = dat[typ]
        return entity_dict,seq_dict

    def angle_order_parameter(self):
        entity_dict, seq_dict = self.get_seq_from_cif()
        for ent in entity_dict.keys():
            if entity_dict[ent] == 'polypeptide(L)':
                ent_index = ent -1
                seq_dat = seq_dict[ent_index]
                for seq in seq_dat.keys():
                    res = seq_dat[seq]
                    angle_list=self.dihedral_angles[res]
                    op_list = []
                    for ang in angle_list:
                        #if seq > min(seq_dat.keys()) and seq < max(seq_dat.keys()):
                        if angle_list.index(ang)==0 and seq > min(seq_dat.keys()):
                            atms = [('{}'.format(seq-1),'{}'.format(ent),seq_dat[seq-1],ang[0]),
                                    ('{}'.format(seq), '{}'.format(ent), res, ang[1]),
                                    ('{}'.format(seq),'{}'.format(ent),res,ang[2]),
                                    ('{}'.format(seq),'{}'.format(ent),res,ang[3])]
                        elif angle_list.index(ang) == 1 and seq < max(seq_dat.keys()):
                            atms = [('{}'.format(seq), '{}'.format(ent), res, ang[0]),
                                    ('{}'.format(seq), '{}'.format(ent), res, ang[1]),
                                    ('{}'.format(seq), '{}'.format(ent), res, ang[2]),
                                    ('{}'.format(seq+1), '{}'.format(ent), seq_dat[seq+1], ang[3])]
                        elif angle_list.index(ang)>1:
                            atms = [('{}'.format(seq), '{}'.format(ent), res, ang[0]),
                                    ('{}'.format(seq), '{}'.format(ent), res, ang[1]),
                                    ('{}'.format(seq), '{}'.format(ent), res, ang[2]),
                                    ('{}'.format(seq), '{}'.format(ent), res, ang[3])]
                        else:
                            atms=None
                        #print (res,atms)
                        if atms is not None:
                            s1 = []
                            for model in range(1,self.max_models+1):
                                coord = self.get_coordinates(model)
                                angl = self.get_dihedral_angle(coord[atms[0]],coord[atms[1]],coord[atms[2]],coord[atms[3]])
                                s1.append(exp(1j*(pi/180)*angl))
                            s = (1.0/len(s1))*abs(sum(s1))
                        else:
                            s = None
                        op_list.append(s)
                    print (seq,res,op_list)







                        #for model in range(1,self.max_models+1)

    @staticmethod
    def get_dihedral_angle(c1, c2, c3, c4):
        # For phi of ith residue provide i-1 C, i N, i CA, i C
        # For psi of ith residue provide i N, i CA, i C, i+1 N

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
        asymid = colnames.index('label_entity_id')
        entityid = colnames.index(
            'label_entity_id')  # asymid is used instead of entity id because of different definition of cif and star
        seqid = colnames.index('label_seq_id')
        co = {}
        for dat in atom_site.getRowList():
            if int(dat[modelid]) == modelID:
                co[(dat[seqid], dat[asymid], dat[compid], dat[atomid])] = numpy.array(
                    [float(dat[xid]), float(dat[yid]), float(dat[zid])])
        return co

if __name__ == '__main__':
    Cyrange('test/data/1nk2.cif')