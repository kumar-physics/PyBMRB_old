#!/usr/bin/env python

from __future__ import print_function

import sys
import ntpath
import os
import json
from numpy import mean,std


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
if PY3:
    from urllib.request import urlopen, Request
else:
    from urllib2 import urlopen, Request
    from cStringIO import StringIO

    BytesIO = StringIO

#http://webapi.bmrb.wisc.edu/v2/list_entries?database=macromolecules
class CS_Prediction(object):
    amino_acids_code_file = scriptPath + '/lib/codeDict.json'

    def __init__(self,atom,nn):
        self.get_entry_list()
        #self.tri_peptide(res, atom)
        #self.tri_peptide_pobabilites(res,atom)
        #self.penta_peptide(atom)
        #self.penta_peptide_pobabilites(atom)
        #self.hepta_peptide(atom)
        #self.hepta_peptide_pobabilites(atom)
        self.random_coil_cs_with_nn(atom,nn)
        self.random_coil_pp(atom,nn)

    @staticmethod
    def get_entry_list():
        url = Request('http://webapi.bmrb.wisc.edu/v2/list_entries?database=macromolecules')
        url.add_header('Application', 'BMRBiViz')
        r = urlopen(url)
        d = json.loads(r.read())
        return d

    @staticmethod
    def get_entry(entryid):
        """
        Downloads the chemical shift data for a given entry id or list of entry ids
        :param entryid: entry or entry ids a list
        :return: chemical shift data
        """

        outdata = []
        if type(entryid) is list:
            for eid in entryid:
                try:
                    indata = pynmrstar.Entry.from_database(eid)
                    cs_data = indata.get_tags(
                        ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                         '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID',
                         '_Atom_chem_shift.Val'])

                    eids = [eid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                    eid_cs_data = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
                                   cs_data['_Atom_chem_shift.Comp_ID'],
                                   cs_data['_Atom_chem_shift.Atom_ID'],
                                   cs_data['_Atom_chem_shift.Atom_type'],
                                   cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                                   cs_data['_Atom_chem_shift.Val']]
                except (OSError, IOError) as e:
                    print(e)
                if len(outdata):
                    for i in range(len(eid_cs_data)):
                        outdata[i] = outdata[i] + eid_cs_data[i]
                else:
                    outdata = eid_cs_data
        else:
            try:
                indata = pynmrstar.Entry.from_database(entryid)
                cs_data = indata.get_tags(
                    ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                     '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID',
                     '_Atom_chem_shift.Val'])
                eids = [entryid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                outdata = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
                           cs_data['_Atom_chem_shift.Comp_ID'],
                           cs_data['_Atom_chem_shift.Atom_ID'],
                           cs_data['_Atom_chem_shift.Atom_type'],
                           cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                           cs_data['_Atom_chem_shift.Val']]
            except (OSError, IOError) as e:
                print(e)
        return outdata

    def tri_peptide(self,res, atom):
        fo = open('tri_peptide_{}_{}.txt'.format(res,atom),'w')
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        aa_dict={'ILE':'I','GLN':'Q','GLY':'G','GLU':'E','CYS':'C','ASP':'D','SER':'S','LYS':'K','PRO':'P','ASN':'N',
                  'VAL':'V','THR':'T','HIS':'H','TRP':'W','PHE':'F','ALA':'A','MET':'M','LEU':'L','ARG':'R','TYR':'Y'}
        tri_pp = []
        for i in aa:
            for j in aa:
                for k in aa:
                    tri_pp.append('{}-{}-{}'.format(j,i,k))

        entry_list = self.get_entry_list()
        for i in entry_list:
            cs_data = self.get_entry(i)
            seq = self.extract_sequence(cs_data)
            #print (sorted(seq.keys()))
            for j in range(len(cs_data[0])):
                if min(seq.keys()) < int(cs_data[1][j]) < max(seq.keys()):
                    if cs_data[2][j] == res and cs_data[3][j] == atom:
                        s = int(cs_data[1][j])
                        try:
                            fo.write('{}-{}-{}\t{}\n'.format(seq[s-1],seq[s],seq[s+1],cs_data[6][j]))
                        except KeyError:
                            pass
    def penta_peptide(self,atom):
        fo = open('penta_peptide_{}.txt'.format(atom),'w')
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        aa_dict={'ILE':'I','GLN':'Q','GLY':'G','GLU':'E','CYS':'C','ASP':'D','SER':'S','LYS':'K','PRO':'P','ASN':'N',
                  'VAL':'V','THR':'T','HIS':'H','TRP':'W','PHE':'F','ALA':'A','MET':'M','LEU':'L','ARG':'R','TYR':'Y'}


        entry_list = self.get_entry_list()
        for i in entry_list:
            cs_data = self.get_entry(i)
            seq = self.extract_sequence(cs_data)
            #print (sorted(seq.keys()))
            for j in range(len(cs_data[0])):
                if min(seq.keys())+1 < int(cs_data[1][j]) < max(seq.keys())-1:
                    if cs_data[3][j] == atom:
                        s = int(cs_data[1][j])
                        try:
                            fo.write('{}-{}-{}-{}-{}\t{}\n'.format(seq[s-2],seq[s-1],seq[s],seq[s+1],seq[s+2],cs_data[6][j]))
                        except KeyError:
                            pass

    def hepta_peptide(self,atom):
        fo = open('hepta_peptide_{}.txt'.format(atom),'w')
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        aa_dict={'ILE':'I','GLN':'Q','GLY':'G','GLU':'E','CYS':'C','ASP':'D','SER':'S','LYS':'K','PRO':'P','ASN':'N',
                  'VAL':'V','THR':'T','HIS':'H','TRP':'W','PHE':'F','ALA':'A','MET':'M','LEU':'L','ARG':'R','TYR':'Y'}


        entry_list = self.get_entry_list()
        for i in entry_list:
            cs_data = self.get_entry(i)
            seq = self.extract_sequence(cs_data)
            #print (sorted(seq.keys()))
            for j in range(len(cs_data[0])):
                if min(seq.keys())+2 < int(cs_data[1][j]) < max(seq.keys())-2:
                    if cs_data[3][j] == atom:
                        s = int(cs_data[1][j])
                        try:
                            fo.write('{}-{}-{}-{}-{}-{}-{}\t{}\n'.format(seq[s-3],seq[s-2],seq[s-1],seq[s],seq[s+1],seq[s+2],seq[s+3],cs_data[6][j]))
                        except KeyError:
                            pass

    def random_coil_cs_with_nn(self,atom,nn=0):
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        fo = open('nn_cs_{}_{}.txt'.format(atom,nn),'w')
        entry_list = self.get_entry_list()
        for entry in entry_list:
            cs_data = self.get_entry(entry)
            seq = self.extract_sequence(cs_data)
            for i in range(len(cs_data[0])):
                if min(seq.keys()) + nn < int(cs_data[1][i]) < max(seq.keys()) - nn:
                    if cs_data[3][i] == atom:
                        s = int(cs_data[1][i])
                        try:
                            if nn == 0 and seq[s] in aa:
                                fo.write('{}\t{}\n'.format(seq[s],cs_data[6][i]))
                            elif nn == 1 and seq[s] in aa and seq[s-1] in aa and seq[s+1] in aa:
                                fo.write('{}-{}-{}\t{}\n'.format(seq[s-1],seq[s],seq[s+1],cs_data[6][i]))
                            elif nn ==2 and seq[s] in aa and seq[s-1] in aa and seq[s+1] in aa and seq[s-2] in aa \
                                    and seq[s+2] in aa:
                                fo.write('{}-{}-{}-{}-{}\t{}\n'.format(seq[s - 2],seq[s - 1], seq[s], seq[s + 1],
                                                                       seq[s + 2],cs_data[6][i]))
                            elif nn == 3 and seq[s] in aa and seq[s-1] in aa and seq[s+1] in aa and seq[s-2] in aa \
                                    and seq[s+2] in aa and seq[s-3] in aa and seq[s+3] in aa:
                                fo.write(
                                    '{}-{}-{}-{}-{}-{}-{}\t{}\n'.format(seq[s - 3],seq[s - 2], seq[s - 1], seq[s],
                                                                        seq[s + 1],seq[s + 2], seq[s + 3],cs_data[6][i]))
                            elif nn == 4 and seq[s] in aa and seq[s-1] in aa and seq[s+1] in aa and seq[s-2] in aa and \
                                    seq[s+2] in aa and seq[s-3] in aa and seq[s+3] in aa and seq[s-4] in aa \
                                    and seq[s+4] in aa:
                                fo.write(
                                    '{}-{}-{}-{}-{}-{}-{}-{}-{}\t{}\n'.format(seq[s - 3],seq[s - 3], seq[s - 2], seq[s - 1],
                                                                              seq[s],seq[s + 1], seq[s + 2], seq[s + 3],
                                                                              seq[s + 4],cs_data[6][i]))
                            else:
                                print ("Not defined !")
                        except KeyError:
                            pass




    def tri_peptide_pobabilites(self,res,atom):
        fo = open('tri_peptide_prediction_{}_{}.txt'.format(res,atom), 'w')
        f = open('tri_peptide_{}_{}.txt'.format(res,atom), 'r').read().split("\n")
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        aa_dict = {'ILE': 'I', 'GLN': 'Q', 'GLY': 'G', 'GLU': 'E', 'CYS': 'C',
                   'ASP': 'D', 'SER': 'S', 'LYS': 'K','PRO': 'P', 'ASN': 'N',
                   'VAL': 'V', 'THR': 'T', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
                   'ALA': 'A', 'MET': 'M', 'LEU': 'L','ARG': 'R', 'TYR': 'Y'}
        tri_pp = {}
        for i in aa:
            for j in aa:
                for k in aa:
                    tri_pp['{}-{}-{}'.format(j, i, k)] = []
        for l in f:
            d=l.split("\t")
            try:
                tri_pp[d[0]].append(float(d[1]))
            except KeyError:
                pass
        for k in tri_pp.keys():
            fo.write('{}\t{}\t{}\n'.format(k,mean(tri_pp[k]),len(tri_pp[k])))
        fo.close()

    def penta_peptide_pobabilites(self,atom):
        fo = open('penta_peptide_prediction_{}.txt'.format(atom), 'w')
        f = open('penta_peptide_{}.txt'.format(atom), 'r').read().split("\n")
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        aa_dict = {'ILE': 'I', 'GLN': 'Q', 'GLY': 'G', 'GLU': 'E', 'CYS': 'C',
                   'ASP': 'D', 'SER': 'S', 'LYS': 'K','PRO': 'P', 'ASN': 'N',
                   'VAL': 'V', 'THR': 'T', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
                   'ALA': 'A', 'MET': 'M', 'LEU': 'L','ARG': 'R', 'TYR': 'Y'}
        tri_pp = {}
        for i in aa:
            for j in aa:
                for k in aa:
                    for l in aa:
                        for m in aa:
                            tri_pp['{}-{}-{}-{}-{}'.format(j, k,i, l,m)] = []
        for l in f:
            d=l.split("\t")
            try:
                tri_pp[d[0]].append(float(d[1]))
            except KeyError:
                pass
        for k in tri_pp.keys():
            fo.write('{}\t{}\t{}\n'.format(k,mean(tri_pp[k]),len(tri_pp[k])))
        fo.close()



    def hepta_peptide_pobabilites(self,atom):
        fo = open('hepta_peptide_prediction_{}.txt'.format(atom), 'w')
        f = open('hepta_peptide_{}.txt'.format(atom), 'r').read().split("\n")[:-1]
        aa = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS', 'ASP', 'SER', 'LYS', 'PRO', 'ASN',
              'VAL', 'THR', 'HIS', 'TRP', 'PHE', 'ALA', 'MET', 'LEU', 'ARG', 'TYR']
        aa_dict = {'ILE': 'I', 'GLN': 'Q', 'GLY': 'G', 'GLU': 'E', 'CYS': 'C',
                   'ASP': 'D', 'SER': 'S', 'LYS': 'K','PRO': 'P', 'ASN': 'N',
                   'VAL': 'V', 'THR': 'T', 'HIS': 'H', 'TRP': 'W', 'PHE': 'F',
                   'ALA': 'A', 'MET': 'M', 'LEU': 'L','ARG': 'R', 'TYR': 'Y'}
        tri_pp = {}
        # for i in aa:
        #     for j in aa:
        #         for k in aa:
        #             for l in aa:
        #                 for m in aa:
        #                     for n in aa:
        #                         for o in aa:
        #                             tri_pp['{}-{}-{}-{}-{}-{}-{}'.format(j,k,l,i,m,n,o)] = []
        for l in f:
            d=l.split("\t")
            #print (d)
            try:
                tri_pp[d[0]].append(float(d[1]))
            except KeyError:
                tri_pp[d[0]] = [float(d[1])]
        for k in tri_pp.keys():
            if len(tri_pp[k]) >0 :
                fo.write('{}\t{}\t{}\n'.format(k,mean(tri_pp[k]),len(tri_pp[k])))
        fo.close()

    def random_coil_pp(self,atom,nn,sd_limt = 8):
        fo = open('nn_pp_{}_{}.txt'.format(atom,nn),'w')
        f = open('nn_cs_{}_{}.txt'.format(atom,nn),'r').read().split("\n")[:-1]
        pep = {}
        for l in f:
            d = l.split("\t")
            try:
                pep[d[0]].append(float(d[1]))
            except KeyError:
                pep[d[0]]= [float(d[1])]
        for k in pep.keys():
            if len(pep[k]) > 0:
                m = mean(pep[k])
                s = std(pep[k])
                min = m - (s*sd_limt)
                max = m + (s*sd_limt)
                filtered = [i for i in pep[k] if min < i < max]
                filt_m = mean(filtered)
                filt_s = std(filtered)
                fo.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(k,m,s,len(pep[k]),filt_m,filt_s,len(filtered)))
        fo.close()


    @staticmethod
    def extract_sequence(cs_data):
        seq = []
        for i in range(len(cs_data[0])):
            seq.append((int(cs_data[1][i]),cs_data[2][i]))
        x=list(set(seq))
        seq_dict = {}
        for i in  x:
            seq_dict[i[0]]=i[1]
        return seq_dict


if __name__ == "__main__":
    atom = sys.argv[1]
    nn = int(sys.argv[2])
    #res = sys.argv[1]
    #atom = sys.argv[2]
    p = CS_Prediction(atom,nn)
