#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import ntpath
import numpy
from math import sqrt,acos
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


class RestraintValidation:

    def __init__(self):
        self.sorted_dist_rest = None
        self.sorted_ang_rest = None
        self.sorted_avg_dist_rest = None
        self.sorted_avg_ang_rest = None
        self.dist_vs_model = []
        self.model_vs_dist = []
        self.ang_vs_model = []
        self.model_vs_ang = []
        self.cal_ang = []
        self.cal_dist = []
        self.ang_rest_info = {}
        self.dist_rest_info = [[], [], [], [], []]
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
        self.dist_rest = None
        self.get_distance_destraints()
        self.dist_restraint_stat()
        print (self.dist_rest_info)
        self.ang_restraint_stat()
        print (self.ang_rest_info)
        self.cal_distance_from_moels()
        self.cal_angle_from_moels()
        self.dist_rest_analysis()
        self.angle_rest_analysis()
        self.sort_angle_violations()
        self.sort_avg_ang_violations()
        self.sort_avg_dist_violations()
        self.sort_dist_violations()


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

    def get_distance_destraints(self):

        """Creates a dictionary of distance restraints from _Gen_dist_constraint loop.
        dist_rest[atom 1,atom 2] = [lower bound, upper bound]
        atoms are identified by (Comp_index_ID,Entity_assembly_ID,Comp_ID,Atom_ID)"""

        dist = self.star_data.get_loops_by_category('_Gen_dist_constraint')  # list of distance restraint loops
        self.dist_rest = []
        for dl in dist:  # for loop in the list
            self.dist_rest.append({})  # each loop info goes into a dictionary
            rest_id = dl.get_tag_names().index("_Gen_dist_constraint.ID")  # rest_id is the key
            atomid1 = dl.get_tag_names().index("_Gen_dist_constraint.Atom_ID_1")  # rest of them are values
            atomid2 = dl.get_tag_names().index("_Gen_dist_constraint.Atom_ID_2")
            compid1 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_ID_1")
            compid2 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_ID_2")
            # entityid1 = dl.get_tag_names().index("_Gen_dist_constraint.Entity_assembly_ID_1")
            # entityid2 = dl.get_tag_names().index("_Gen_dist_constraint.Entity_assembly_ID_2")
            entityid1 = dl.get_tag_names().index("_Gen_dist_constraint.Auth_asym_ID_1")
            # asymid is used instead of entity id because of different definition of cif and star
            entityid2 = dl.get_tag_names().index("_Gen_dist_constraint.Auth_asym_ID_2")
            seqid1 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_index_ID_1")
            seqid2 = dl.get_tag_names().index("_Gen_dist_constraint.Comp_index_ID_2")
            lb = dl.get_tag_names().index("_Gen_dist_constraint.Distance_lower_bound_val")
            ub = dl.get_tag_names().index("_Gen_dist_constraint.Distance_upper_bound_val")
            for i in dl:  # for every row in the restraint loop
                if int(i[rest_id]) not in self.dist_rest[-1].keys():
                    # each rest_id may have multiple rows, because of ambiguity codes and NEF atom
                    # nomeclature ; so it will be a list
                    self.dist_rest[-1][int(i[rest_id])] = []
                try:
                    lbv = float(i[lb])
                except ValueError:
                    lbv = -999
                try:
                    ubv = float(i[ub])
                except ValueError:
                    ubv = 999
                if lbv == -999 and ubv == 999:
                    print (
                        "Distance restraint value not readable {},{} between ({}{}{}{})-({}{}{}{})".format(i[lb], i[ub],
                                                                                                           i[seqid1],
                                                                                                           i[entityid1],
                                                                                                           i[compid1],
                                                                                                           i[atomid1],
                                                                                                           i[seqid2],
                                                                                                           i[entityid2],
                                                                                                           i[compid2],
                                                                                                           i[atomid2]))
                else:
                    self.dist_rest[-1][int(i[rest_id])].append([(i[seqid1], i[entityid1], i[compid1], i[atomid1]),
                                                                (i[seqid2], i[entityid2], i[compid2], i[atomid2]),
                                                                (lbv, ubv)])
                # self.dist_rest[(i[seqid1],i[entityid1],i[compid1],i[atomid1]),(i[seqid2],i[entityid2],i[compid2],
                # i[atomid2])]= [lbv,ubv]
        print ('Number of distance restraints : {}'.format([len(i) for i in self.dist_rest]))

    def dist_restraint_stat(self):
        for rlistid in range(len(self.dist_rest)):
            rest_list = self.dist_rest[rlistid]
            for restid in rest_list.keys():
                for rest in rest_list[restid]:
                    try:
                        n = abs(int(rest[0][0]) - int(rest[1][0]))
                    except ValueError:
                        self.logger.warning("Invalid restraint {}".format(rest))
                self.dist_rest_info[0].append((restid, rlistid))
                if n == 0:
                    self.dist_rest_info[1].append((restid, rlistid))
                elif n == 1:
                    self.dist_rest_info[2].append((restid, rlistid))
                elif n > 1 and n < 5:
                    self.dist_rest_info[3].append((restid, rlistid))
                else:
                    self.dist_rest_info[4].append((restid, rlistid))
        # self.logger.debug('Distance restraint statistics {}'.format(self.dist_rest_info))
        return self.dist_rest_info

    def ang_restraint_stat(self):
        for rest_list in self.ang_rest:
            for rest in rest_list.keys():
                if rest_list[rest][0] not in self.ang_rest_info.keys():
                    self.ang_rest_info[rest_list[rest][0]] = 0
                self.ang_rest_info[rest_list[rest][0]] += 1
        return self.ang_rest_info
        #self.logger.debug('Angle restraint statistics {}'.format(self.ang_rest_info))

    def cal_distance_from_moels(self):
        c = []
        for modelID in range(1, self.max_models + 1):
            c.append(self.get_coordinates(modelID))
        for rest_list in self.dist_rest:
            self.cal_dist.append({})
            for rest_id in rest_list.keys():
                if rest_id not in self.cal_dist[-1].keys():
                    self.cal_dist[-1][rest_id] = []
                for rest in rest_list[rest_id]:
                    c_dist = [rest[2]]
                    for m in c:
                        try:
                            d = self.get_distance(m[rest[0]], m[rest[1]])
                        except KeyError:
                            print (
                                'One of the atom from the restraint not found in coordinates {},{},{}'.format(rest_id,
                                                                                                              rest[0],
                                                                                                              rest[1]))
                            d = '.'
                        c_dist.append(d)
                    self.cal_dist[-1][rest_id].append(c_dist)

    def cal_angle_from_moels(self):
        c = []
        for modelID in range(1, self.max_models + 1):
            c.append(self.get_coordinates(modelID))
        for rest_list in self.ang_rest:
            self.cal_ang.append({})
            for rest_id in rest_list.keys():
                rest = rest_list[rest_id]

                if rest_id not in self.cal_ang[-1].keys():
                    self.cal_ang[-1][rest_id] = [rest[5]]

                for m in c:
                    try:
                        d = self.get_dihedral_angle(m[rest[1]], m[rest[2]], m[rest[3]], m[rest[4]])
                    except KeyError:
                        print (
                            'One of the atom from the restraint not found in coordinates {},{},{},{},{}'.format(rest_id,
                                                                                                                rest[1],
                                                                                                                rest[2],
                                                                                                                rest[2],
                                                                                                                rest[
                                                                                                                    2]))
                        d = '.'
                    self.cal_ang[-1][rest_id].append(d)

    def dist_rest_analysis(self):
        for m in range(self.max_models):
            self.model_vs_dist.append([0, 0, 0, 0, 0, 0, 0, 0])
            for rest_list in self.cal_dist:
                for key in sorted(rest_list):
                    self.model_vs_dist[-1][0] += 1
                    dav = []
                    d_flag = False
                    d_list = []
                    for rest in rest_list[key]:
                        d_list.append(rest[m + 1])

                    r6av = self.r6sum(d_list)
                    if r6av < rest[0][0] or r6av > rest[0][1]:
                        d_flag = True
                        if r6av < rest[0][0]:
                            dav.append(rest[0][0] - r6av)
                        else:
                            dav.append(r6av - rest[0][1])

                    #                         d_min = rest[0][0]
                    #                         d_max = rest[0][1]
                    #                         d=rest[m+1]
                    #                         if d>d_max or d<d_min:
                    #                             d_flag=True
                    #                             if d<d_max:
                    #                                 dd=d_max-d
                    #                             else:
                    #                                 dd=d-d_min
                    #                             dav.append(dd)
                    if d_flag:
                        dif = sum(dav) / len(dav)

                        if dif > 0 and dif <= 0.2:
                            self.model_vs_dist[-1][1] += 1
                        elif dif > 0.2 and dif <= 0.5:
                            self.model_vs_dist[-1][2] += 1
                        elif dif > 0.5 and dif <= 1.0:
                            self.model_vs_dist[-1][3] += 1
                        elif dif > 1.0 and dif <= 2.0:
                            self.model_vs_dist[-1][4] += 1
                        elif dif > 2.0 and dif <= 5.0:
                            self.model_vs_dist[-1][5] += 1
                        else:
                            self.model_vs_dist[-1][6] += 1
                    else:
                        self.model_vs_dist[-1][7] += 1
        for rest_list in self.cal_dist:
            self.dist_vs_model.append([])
            for key in sorted(rest_list):
                mm = [key]
                # print rest
                for i in range(1, self.max_models + 1):
                    d_flag = False
                    err = []
                    d_list = []
                    for rest in rest_list[key]:
                        d_list.append(rest[i])
                    r6avg = self.r6sum(d_list)
                    if r6avg < rest[0][0] or r6avg > rest[0][1]:
                        d_flag = True
                        if r6avg < rest[0][0]:
                            err.append(rest[0][0] - r6avg)
                        else:
                            err.append(r6avg - rest[0][1])
                    #                         if rest[i]<rest[0][0] or rest[i]>rest[0][1]:
                    #                             d_flag=True
                    #                             if rest[i]<rest[0][0]:
                    #                                 err.append(rest[0][0]-rest[i])
                    #                             else:
                    #                                 err.append(rest[i]-rest[0][1])
                    if d_flag:
                        err_avg = sum(err) / len(err)
                        mm.append((i, err_avg))

                self.dist_vs_model[-1].append(mm)
        return [self.model_vs_dist, self.dist_vs_model]

    def angle_rest_analysis(self):
        for m in range(self.max_models):
            self.model_vs_ang.append([0, 0, 0, 0, 0, 0, 0, 0])
            for rest_list in self.cal_ang:
                for key in sorted(rest_list):
                    self.model_vs_ang[-1][0] += 1
                    dav = []
                    d_flag = False
                    rest = rest_list[key]
                    d_min = rest[0][0]
                    d_max = rest[0][1]
                    d = rest[m + 1]
                    dif = 0.0
                    if (d_max > 0 and d_min > 0) and d < 0:
                        d = 360.0 + d
                    if (d_max < 0 and d_min < 0) and d > 0:
                        d = -1 * (360.0 - d)
                    if d > d_max:
                        dif = d - d_max
                        d_flag = True
                        #print ('dmax {}:{},{},{},{},{}'.format(key, d_min, d_max, d, dif, d_flag))
                    elif d < d_min:
                        dif = d_min - d
                        d_flag = True
                        #print ('dmin {}:{},{},{},{},{}'.format(key, d_min, d_max, d, dif, d_flag))
                    else:
                        d_flag = False
                        #print ('dmin {}:{},{},{},{},{}'.format(key, d_min, d_max, d, dif, d_flag))
                    #print ('{}:{},{},{},{}'.format(key, d_min, d_max, d, dif))
                    if d_flag:
                        if dif > 0 and dif <= 5:
                            self.model_vs_ang[-1][1] += 1
                        elif dif > 5 and dif <= 10:
                            self.model_vs_ang[-1][2] += 1
                        elif dif > 10 and dif <= 20:
                            self.model_vs_ang[-1][3] += 1
                        elif dif > 20 and dif <= 40:
                            self.model_vs_ang[-1][4] += 1
                        elif dif > 40 and dif <= 80:
                            self.model_vs_ang[-1][5] += 1
                        else:
                            self.model_vs_ang[-1][6] += 1
                    else:
                        self.model_vs_ang[-1][7] += 1
        for rest_list in self.cal_ang:
            self.ang_vs_model.append([])
            for key in sorted(rest_list):
                mm = [key]
                # print rest
                for i in range(1, self.max_models + 1):
                    d_flag = False
                    rest = rest_list[key]
                    rest_ang = rest[i]
                    if (rest[0][0] > 0 and rest[0][1] > 0) and rest[i] < 0:
                        rest_ang = 360.0 + rest[i]
                    if (rest[0][0] < 0 and rest[0][1] < 0) and rest[i] > 0:
                        rest_ang = -1 * (360.0 - rest[i])
                    if rest_ang < rest[0][0] or rest_ang > rest[0][1]:
                        d_flag = True

                        if rest_ang < rest[0][0]:
                            err_ang = rest[0][0] - rest_ang
                        else:
                            err_ang = rest_ang - rest[0][1]
                    if d_flag:
                        mm.append((i, err_ang))
                self.ang_vs_model[-1].append(mm)
        return [self.model_vs_ang, self.ang_vs_model]

    def sort_dist_violations(self):
        r_list = []
        for j in self.dist_vs_model:
            for i in j:
                if len(i) > 1:
                    for k in range(1, len(i)):
                        r_list.append([i[0], i[k][0], i[k][1]])
        self.sorted_dist_rest = sorted(r_list, key=operator.itemgetter(2), reverse=True)

    def sort_avg_dist_violations(self):
        r_list = []
        for j in self.dist_vs_model:
            for i in j:
                if len(i) > 1:
                    rid = i[0]
                    err = [x[1] for x in i[1:]]
                    r_list.append([rid, len(err), sum(err) / len(err)])

        self.sorted_avg_dist_rest = sorted(r_list, key=operator.itemgetter(2), reverse=True)

    def sort_avg_ang_violations(self):
        r_list = []
        for j in self.ang_vs_model:
            for i in j:
                if len(i) > 1:
                    rid = i[0]
                    err = [x[1] for x in i[1:]]
                    r_list.append([rid, len(err), sum(err) / len(err)])

        self.sorted_avg_ang_rest = sorted(r_list, key=operator.itemgetter(2), reverse=True)

    def sort_angle_violations(self):
        r_list = []

        for j in self.ang_vs_model:
            for i in j:
                if len(i) > 1:
                    for k in range(1, len(i)):
                        r_list.append([i[0], i[k][0], i[k][1]])
        self.sorted_ang_rest = sorted(r_list, key=operator.itemgetter(2), reverse=True)

if __name__ == "__main__":
    p = RestraintValidation()
