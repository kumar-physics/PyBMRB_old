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
        self.centroid_shifted, self.mean_coord = self.cal_mean_coord(['N', 'CA', 'C'])
        self.cluster_atoms()
        #print (self.cal_rmsd(self.core_atoms,0))


    def cal_mean_coord(self,atm_list):
        x = []
        c=[]
        for model in range(self.max_models):
            y=[]
            c1=[]
            for k in self.coord_data[model].keys():
                if k[3] in atm_list:
                    c1.append(k)
                    y.append(self.coord_data[model][k])
            x.append(numpy.array(y))
            c.append(c1)
        m = numpy.empty(numpy.shape(x[0]))

        for p in x:
            pc = self.centroid(p)
            p-=pc

        for p in x:
            P1 = deepcopy(p)
            #pc = self.centroid(P1)
            #P1 -= pc
            m += P1
        m = (1.0 / len(x)) * m
        mean_coord = {}
        for i in range(len(m)):
            mean_coord[c[0][i]] = m[i]
        centroid_coord = []
        for i in range(len(x)):
            d = {}
            for j in range(len(x[i])):
                d[c[i][j]] = x[i][j]
            centroid_coord.append(d)
        return centroid_coord, mean_coord







    def cal_rmsd(self,atmlist,ent_id):
        #print (self.seq_dict)
        alist = []
        for atm in atmlist:
            #print (atm)
            alist.append(('{}'.format(atm),'{}'.format(ent_id+1),'{}'.format(self.seq_dict[ent_id][int(atm)]),'N'))
            alist.append(('{}'.format(atm), '{}'.format(ent_id+1), '{}'.format(self.seq_dict[ent_id][int(atm)]), 'CA'))
            alist.append(('{}'.format(atm), '{}'.format(ent_id+1), '{}'.format(self.seq_dict[ent_id][int(atm)]), 'C'))

        x=[]
        for model in range(self.max_models):
            y=[]
            for atm in alist:
                y.append(self.centroid_shifted[model][atm])
            x.append(numpy.array(y))
        m = []
        for atm in alist:
            m.append(self.mean_coord[atm])


        # m = numpy.empty(numpy.shape(x[0]))
        # print (atmlist,len(m))
        # for p in x:
        #     P1 = deepcopy(p)
        #     pc = self.centroid(P1)
        #     P1-=pc
        #     m+=P1
        # m=(1.0/len(x))*m
        r = 0.0
        for p in x:
            #P1 = deepcopy(p)
            #pc = self.centroid(P1)
            #P1 -= pc
            r+=self.kabsch(m,p)
        return (1.0/(len(x)))*r


    def calculate_a(self,cluster_list,ent_id):
        r=0.0
        n=0
        for c in cluster_list:
            n+=len(c)
            r += self.cal_rmsd(c,ent_id)
        return (1.0/n)*r





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
                    s = 0.0
                    for k in range(len(distance_matrix)):
                        s += (distance_matrix[k][(i, j)] - mean_distance[(i, j)]) ** 2.0
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

    @staticmethod
    def order_seq(seqlist):
        seqint = [int(i) for i in seqlist]
        return ['{}'.format(i) for i in sorted(seqint)]

    def cluster_atoms(self):
        for ent in self.d_matrix:
            if ent is not None:
                dmat = ent
                atm_list = [i[0] for i in self.core_atoms[self.d_matrix.index(ent)]]
                cluster_list=[]
                n = len(atm_list)
                a=[]
                ni = []
                for i in range(n):
                    #print (len(cluster_list),[len(o) for o in cluster_list])
                    if len(cluster_list) == 0:
                        v = self.v_value(atm_list,dmat)
                        vmin_idx = min(v, key=v.get)
                        vmin_value = v[vmin_idx]
                        cluster_list.append([j for j in vmin_idx])
                        #print (vmin_idx,vmin_value,cluster_list)
                        atm_list = [j for j in atm_list if j not in vmin_idx]
                    else:
                        a.append(self.calculate_a(cluster_list,self.d_matrix.index(ent)))
                        ni.append(len(cluster_list)+len(atm_list))
                        if len(atm_list) > 1:
                            v = self.v_value(atm_list, dmat)
                            vmin_idx = min(v, key=v.get)
                            vmin_value = v[vmin_idx]
                        else:
                            vmin_idx =('X')
                            vmin_value = 9999999.999
                        vc = {}
                        if len(atm_list) > 0:
                            for atm in atm_list:
                                for c in cluster_list:
                                    x=deepcopy(c)
                                    x.append(atm)
                                    v2 = self.v_value(x,dmat)
                                    vv=[v2[k] for k in v2.keys()]
                                    vc[tuple(x)] = variance(vv)
                            vcmin_idx = min(vc, key=vc.get)
                            vcmin_value = vc[vcmin_idx]
                        else:
                            vcmin_idx =('X')
                            vcmin_value = 9999999.999
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
                            vvcmin_value = 9999999.999
                        if vmin_value < vcmin_value and vmin_value < vvcmin_value:
                            cluster_list.append([j for j in self.order_seq(vmin_idx)])
                            atm_list = [j for j in atm_list if j not in vmin_idx]
                        elif vcmin_value < vmin_value and vcmin_value < vvcmin_value:
                            for k in vcmin_idx:
                                for k2 in cluster_list:
                                    if k in k2:
                                        cluster_list.remove(k2)
                            cluster_list.append([j for j in self.order_seq(vcmin_idx)])
                            atm_list = [j for j in atm_list if j not in vcmin_idx]
                        elif vvcmin_value < vmin_value and vvcmin_value < vcmin_value:
                            for k in vvcmin_idx:
                                for k2 in cluster_list:
                                    if k in k2:
                                        cluster_list.remove(k2)
                            cluster_list.append([j for j in self.order_seq(vvcmin_idx)])
                            atm_list = [j for j in atm_list if j not in vvcmin_idx]
                print (n,len(a),len(ni))
                for i in range(n-1):
                    amin = min(a)
                    amax = max(a)
                    p = ((n-1)*(a[i]-amin)/(amax-amin))+ni[i]
                    print (1,p)










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
            if int(dat[modelid]) == modelID: #and dat[asymid]=='polypeptide(L)':
                co[(dat[seqid], dat[asymid], dat[compid], dat[atomid])] = numpy.array(
                    [float(dat[xid]), float(dat[yid]), float(dat[zid])])
        return co

    def fit(self,P, Q):
        """ Varies the distance between P and Q, and optimizes rotation for each step
        until a minimum is found.
        """
        step_size = P.max(0)
        threshold = step_size * 1e-9
        rmsd_best = self.kabsch(P, Q)
        while True:
            for i in range(3):
                temp = numpy.zeros(3)
                temp[i] = step_size[i]
                rmsd_new = self.kabsch(P + temp, Q)
                if rmsd_new < rmsd_best:
                    rmsd_best = rmsd_new
                    P[:, i] += step_size[i]
                else:
                    rmsd_new = self.kabsch(P - temp, Q)
                    if rmsd_new < rmsd_best:
                        rmsd_best = rmsd_new
                        P[:, i] -= step_size[i]
                    else:
                        step_size[i] /= 2
            if (step_size < threshold).all():
                break
        return rmsd_best

    def kabsch(self,P, Q):
        """ The Kabsch algorithm

        http://en.wikipedia.org/wiki/Kabsch_algorithm

        The algorithm starts with two sets of paired points P and Q.
        P and Q should already be centered on top of each other.

        Each vector set is represented as an NxD matrix, where D is the
        the dimension of the space.

        The algorithm works in three steps:
        - a translation of P and Q
        - the computation of a covariance matrix C
        - computation of the optimal rotation matrix U

        The optimal rotation matrix U is then used to
        rotate P unto Q so the RMSD can be caculated
        from a straight forward fashion.

        """

        # Computation of the covariance matrix
        C = numpy.dot(numpy.transpose(P), Q)

        # Computation of the optimal rotation matrix
        # This can be done using singular value decomposition (SVD)
        # Getting the sign of the det(V)*(W) to decide
        # whether we need to correct our rotation matrix to ensure a
        # right-handed coordinate system.
        # And finally calculating the optimal rotation matrix U
        # see http://en.wikipedia.org/wiki/Kabsch_algorithm
        V, S, W = numpy.linalg.svd(C)
        d = (numpy.linalg.det(V) * numpy.linalg.det(W)) < 0.0

        if (d):
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # Create Rotation matrix U
        U = numpy.dot(V, W)

        # Rotate P
        P = numpy.dot(P, U)
        return self.rmsd(P, Q)
    @staticmethod
    def centroid(X):
        """ Calculate the centroid from a vectorset X """
        C = sum(X) / len(X)
        return C
    @staticmethod
    def rmsd(V, W):
        """ Calculate Root-mean-square deviation from two sets of vectors V and W.
        """
        D = len(V[0])
        N = len(V)
        rmsd = 0.0
        for v, w in zip(V, W):
            rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
        return numpy.sqrt(rmsd / N)

if __name__ == '__main__':
    Cyrange('test/data/2mtv.cif')