#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import ntpath
import numpy
from math import sqrt,acos,pi
from cmath import exp
from copy import deepcopy
from statistics import variance

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
        self.cif_data, self.max_models, self.entity_dict, self.seq_dict = self.read_ciffile(ciffile)
        self.coord_data = self.get_models()
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
        self.op_dict = self.angle_order_parameter()
        self.scutoff = self.calculate_scutoff(self.op_dict)
        self.core_atoms = self.calculate_core_atoms(self.op_dict,self.scutoff,self.seq_dict)
        self.d_matrix=self.distance_matrix()
        self.cluster_atoms()

    def get_models(self):
        coord = []
        for model in range(1,self.max_models+1):
            coord.append(self.get_coordinates(model))
        return coord

    @staticmethod
    def read_ciffile(ciffile):
        """ Reads the input coordinate CIF file using pdbx lightweight parser
        ref: http://mmcif.wwpdb.org/docs/sw-examples/python/html/index.html """
        cif_data = []
        ifh = open(ciffile, 'r')
        pRd = PdbxReader(ifh)
        pRd.read(cif_data)
        ifh.close()
        c0 = cif_data[0]
        atom_site = c0.getObj('atom_site')
        max_models = int(atom_site.getValue('pdbx_PDB_model_num', -1))
        seq_dat = c0.getObj('entity_poly_seq')
        colnames = seq_dat.getAttributeList()
        entity_id = colnames.index('entity_id')
        seq_id = colnames.index('num')
        res_id = colnames.index('mon_id')
        max_entity = int(seq_dat.getValue('entity_id', -1))
        seq_dict = []
        for i in range(max_entity):
            seq_dict.append({})
        for dat in seq_dat.getRowList():
            seq_dict[int(dat[entity_id]) - 1][int(dat[seq_id])] = dat[res_id]
        ent_dat = c0.getObj('entity_poly')
        colnames = ent_dat.getAttributeList()
        entity_id = colnames.index('entity_id')
        typ = colnames.index('type')
        entity_dict = {}
        for dat in ent_dat.getRowList():
            entity_dict[int(dat[entity_id])] = dat[typ]
        return cif_data,max_models,entity_dict,seq_dict

    # def get_seq_from_cif(self):
    #     self.cif_data = []
    #     ifh = open(self.cfile, 'r')
    #     pRd = PdbxReader(ifh)
    #     pRd.read(self.cif_data)
    #     ifh.close()
    #     c0 = self.cif_data[0]
    #     seq_dat = c0.getObj('entity_poly_seq')
    #     colnames = seq_dat.getAttributeList()
    #     entity_id = colnames.index('entity_id')
    #     seq_id = colnames.index('num')
    #     res_id = colnames.index('mon_id')
    #     max_entity = int(seq_dat.getValue('entity_id',-1))
    #     seq_dict=[]
    #     for i in range(max_entity):
    #         seq_dict.append({})
    #     for dat in seq_dat.getRowList():
    #         seq_dict[int(dat[entity_id])-1][int(dat[seq_id])]=dat[res_id]
    #     ent_dat = c0.getObj('entity_poly')
    #     colnames = ent_dat.getAttributeList()
    #     entity_id = colnames.index('entity_id')
    #     typ = colnames.index('type')
    #     entity_dict = {}
    #     for dat in ent_dat.getRowList():
    #         entity_dict[int(dat[entity_id])] = dat[typ]
    #     return entity_dict,seq_dict

    def angle_order_parameter(self):
        entity_dict= self.entity_dict
        seq_dict = self.seq_dict
        op_dict=[]
        for ent in entity_dict.keys():
            op_dict.append({})
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
                            for model in range(self.max_models):
                                coord = self.coord_data[model]
                                angl = self.get_dihedral_angle(coord[atms[0]],coord[atms[1]],coord[atms[2]],coord[atms[3]])
                                s1.append(exp(1j*(pi/180)*angl))
                            s = (1.0/len(s1))*abs(sum(s1))
                        else:
                            s = None
                        op_list.append(s)
                    op_dict[ent_index][seq] = op_list
        return op_dict

    @staticmethod
    def calculate_scutoff(op_dict):
        scutoff = []
        for ent in op_dict:
            if len(ent):
                svalues = []
                for k in ent.keys():
                    svalues+=ent[k]
                svalues = [i for i in svalues if i is not None]
                sorted_svalues = sorted(svalues)
                q=[]
                for i in range(len(sorted_svalues)):
                    q.append(((len(sorted_svalues) - 1) * ((sorted_svalues[i] - min(sorted_svalues)) /
                                                           (max(sorted_svalues) - min(sorted_svalues)))) - (i + 1))
                qindex = q.index(max(q))
                scutoff.append(sorted_svalues[qindex])
            else:
                scutoff.append(None)
        return scutoff

    @staticmethod
    def calculate_core_atoms(op_dict, scutoff, seq_dict):
        core_atoms = []
        for ent in op_dict:
            if len(ent):
                core_atoms.append([])
                for k in ent.keys():
                    check_cutoff = [i > scutoff[op_dict.index(ent)] for i in ent[k] if i is not None]
                    if True in check_cutoff:
                        catom = ('{}'.format(k), '{}'.format(op_dict.index(ent) + 1), seq_dict[op_dict.index(ent)][k], 'CA')
                        core_atoms[op_dict.index(ent)].append(catom)
            else:
                core_atoms.append(None)
        return core_atoms

    def distance_matrix(self):
        v=[]
        for ent in self.core_atoms:
            if ent is not None:
                y=[]
                for model in range(self.max_models):
                    coord = self.coord_data[model]
                    m={}
                    for atm1 in ent:
                        for atm2 in ent:
                            if atm1 != atm2:
                                dist = self.get_distance(coord[atm1],coord[atm2])
                                m[(atm1[0],atm2[0])]=dist
                    y.append(m)
                v.append(y)
            else:
                v.append(None)
        return v

    def v_value(self,atom_list, distance_matrix,clust_list=None):
        mean_distance = {}
        for i in atom_list:
            for j in atom_list:
                if i != j:
                    s = 0.0
                    for k in range(len(distance_matrix)):
                        s += distance_matrix[k][(i, j)]
                    mean_distance[(i, j)] = (1.0 / len(distance_matrix)) * s
        v = {}
        for i in atom_list:
            for j in atom_list:
                if i != j and self.not_within_cluster(i,j,clust_list):
                    s = 0
                    for k in range(len(distance_matrix)):
                        s += (distance_matrix[k][(i, j)] - mean_distance[(i, j)]) ** 2
                    v[(i, j)] = (1.0 / len(distance_matrix)) * s
        return v

    @staticmethod
    def not_within_cluster(i, j, cluster_list):
        if cluster_list is None:
            flg = True
        else:
            flg = True
            for c in cluster_list:
                if i in c and j in c:
                    flg = False
        return flg

    def cluster_atoms(self):
        for ent in self.d_matrix:
            if ent is not None:
                dmat = ent
                atm_list = [i[0] for i in self.core_atoms[self.d_matrix.index(ent)]]
                cluster_list=[]
                n = len(atm_list)
                for i in range(n):
                    print (cluster_list)
                    if len(cluster_list) == 0:
                        v = self.v_value(atm_list,dmat)
                        vmin_idx = min(v, key=v.get)
                        vmin_value = v[vmin_idx]
                        cluster_list.append([j for j in vmin_idx])
                        print (vmin_idx,vmin_value,cluster_list)
                        atm_list = [j for j in atm_list if j not in vmin_idx]
                    else:
                        v = self.v_value(atm_list, dmat)
                        vmin_idx = min(v, key=v.get)
                        vmin_value = v[vmin_idx]
                        vc = {}
                        for atm in atm_list:
                            for c in cluster_list:
                                x=deepcopy(c)
                                x.append(atm)
                                v2 = self.v_value(x,dmat)
                                vv=[v2[k] for k in v2.keys()]
                                vc[tuple(x)] = variance(vv)
                        vcmin_idx = min(vc, key=vc.get)
                        vcmin_value = vc[vcmin_idx]
                        vvc = {}
                        if len(cluster_list)>1:
                            for c1 in cluster_list:
                                for c2 in cluster_list:
                                    if c1 != c2:
                                        x = c1+c2
                                        v2 = self.v_value(x, dmat)
                                        vv = [v2[k] for k in v2.keys()]
                                        vvc[tuple(x)] = variance(vv)
                            vvcmin_idx = min(vvc, key=vvc.get)
                            vvcmin_value = vvc[vvcmin_idx]
                        else:
                            vvcmin_idx = ('x')
                            vvcmin_value = 9999.999
                        if vmin_value < vcmin_value and vmin_value < vvcmin_value:
                            cluster_list.append([j for j in vmin_idx])
                            atm_list = [j for j in atm_list if j not in vmin_idx]
                        elif vcmin_value < vmin_value and vcmin_value < vvcmin_value:
                            for k in vcmin_idx:
                                for k2 in cluster_list:
                                    if k in k2:
                                        cluster_list.remove(k2)
                            cluster_list.append([j for j in vcmin_idx])
                            atm_list = [j for j in atm_list if j not in vcmin_idx]
                        elif vvcmin_value < vmin_value and vvcmin_value < vcmin_value:
                            for k in vvcmin_idx:
                                for k2 in cluster_list:
                                    if k in k2:
                                        cluster_list.remove(k2)
                            cluster_list.append([j for j in vvcmin_idx])
                            atm_list = [j for j in atm_list if j not in vvcmin_idx]










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

    @staticmethod
    def get_distance(c1, c2):

        """ Calculates the distance between two coordinate values.
        Each coordinate is an array of x,y,z. Returns distance in A assuming the input coordinates are in A"""

        return numpy.linalg.norm(c1 - c2)

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
    Cyrange('test/data/1cfc.cif')