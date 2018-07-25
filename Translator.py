#!/usr/bin/env python

from __future__ import print_function

import csv
import json
import re
import sys
import os
import ntpath
import datetime
import time


PY3 = (sys.version_info[0] == 3)

(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))

try:
    import pynmrstar

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

__version__ = "v1.0"


class Translator:
    """
    NEF to NMR-STAR translator
    """
    map_file = scriptPath + '/lib/NEF_NMRSTAR_equivalence.csv'
    nef_info_file = scriptPath + '/lib/NEF_mandatory.csv'
    atom_nomenclature_file = scriptPath + '/lib/atomDict.json'
    amino_acids_code_file = scriptPath + '/lib/codeDict.json'

    def __init__(self):
        self.atom_dict = None
        self.code_dict = None
        self.tag_map = None
        self.nef_info = None
        if not self.load_atom_dict():
            err_msg = "Can't create atom nomenclature dictionary." \
                      "File {} either missing or not readable".format(self.atom_nomenclature_file)
            print (err_msg)
        if not self.load_code_dict():
            err_msg = "Can't create amio acids three letter to single letter dictionary. " \
                      "File {} either missing or not readable".format(self.amino_acids_code_file)
            print (err_msg)
        if not self.load_map_file():
            err_msg = "Can't create NEF to NMR-STAR mapping dictionary. " \
                      "File {} either missing or not readable".format(self.map_file)
            print(err_msg)
        if not self.load_nef_info():
            err_msg = "Can't create NEF mandatory tag list. " \
                      "File {} either missing or not readable".format(self.nef_info_file)
            print(err_msg)
        self.chains = None

    @staticmethod
    def read_input(in_file):
        """Reads the input file.Returns file readable(T/F),
        content type(Entry/Saveframe/loop), data
        :param in_file: Input NMR-STAR/NEF file name
        :return: file readble, content type, data (STAR object)
        """
        file_readable = False
        content_type = None
        try:
            in_data = pynmrstar.Entry.from_file(in_file)
            file_readable = True
            content_type = "Entry"
        except ValueError:
            try:
                in_data = pynmrstar.Saveframe.from_file(in_file)
                file_readable = True
                content_type = "Saveframe"
            except ValueError:
                try:
                    in_data = pynmrstar.Loop.from_file(in_file)
                    file_readable = True
                    content_type = "Loop"
                except ValueError as err:
                    in_data = None
                    err_msg = "File contains no valid saveframe or loop. Invalid NMRSTAR file :{}".format(err)
                    print(err_msg)
        except IOError as err:
            in_data = None
            err_msg = "File not found:{}".format(err)
            print(err_msg)
        return file_readable, content_type, in_data

    def load_atom_dict(self):
        """Reads the atomDict.json file and creates a dictionary of residues and atoms"""
        status = False
        try:
            with open(self.atom_nomenclature_file, 'r') as atom_file:
                self.atom_dict = json.loads(atom_file.read())
            status = True
        except IOError:
            err_msg = "atomDict.json file is missing! Check whether the file is inside {}".format(scriptPath + '/lib')
            print(err_msg)
        return status

    def load_code_dict(self):
        status = False
        """Reads the codeDict.json file and creates a dictionary of residues and atoms"""
        try:
            with open(self.amino_acids_code_file, 'r') as code_file:
                self.code_dict = json.loads(code_file.read())
            status = True
        except IOError:
            err_msg = "codeDict.json file is missing! Check whether the file is inside {}".format(scriptPath + '/lib')
            print(err_msg)
        return status

    def get_one_letter_code(self, res):
        """
        Takes 3 letter code as input and return one letter code. If one letter code is not defined in the dictionary
        then it will return ?
        :param res: 3 letter amino acid code
        :type res: str
        :return: single letter amino acid code
        :rtype : str
        """
        try:
            one_letter_code = self.code_dict[res.upper()]
        except KeyError:
            one_letter_code = '?'
        return one_letter_code

    def load_map_file(self):
        """
        Reads the NEF_NMRSTAR_equivalence.csv file and creates mapping as a list
        """
        status = False
        try:
            with open(self.map_file, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                map_data = []
                for r in spamreader:
                    if r[0][0] != '#':
                        map_data.append(r)
            self.tag_map = list(map(list, zip(*map_data)))
            status = True
        except IOError:
            err_msg = "NEF-NMRSTAR_equivalence.csv file is missing! check the file is inside  {} ".format(
                scriptPath + '/lib')
            print(err_msg)
            self.tag_map = []
        return status

    def load_nef_info(self):
        """
        Reads mandatory tag information for NEF file
        """
        status = False
        try:
            with open(self.nef_info_file, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                map_data = []
                for r in spamreader:
                    if r[0][0] != '#':
                        map_data.append(r)
            self.nef_info = map_data
            status = True
        except IOError:
            err_msg = "NEF_mandatory.csv file is missing! check the file is inside  {} ".format(scriptPath + '/lib')
            print(err_msg)
            self.nef_info = []
        return status

    def validate_file(self, input_file):
        valid_cs_file = False
        valid_re_file = False
        file_type = None
        readable, content, input_data = self.read_input(input_file)
        if readable:
            sf_list, lp_list = self.get_saveframes_and_loops(input_data, content)
            if '_nef_chemical_shift' in lp_list and \
                    not self.is_loop_empty(input_data, '_nef_chemical_shift', content):
                valid_cs_file = True
                file_type = "NEF"
            if '_nef_distance_restraint' in lp_list and \
                    not self.is_loop_empty(input_data, '_nef_distance_restraint', content):
                valid_re_file = True
                file_type = "NEF"
            if '_Atom_chem_shift' in lp_list and \
                    not self.is_loop_empty(input_data, '_Atom_chem_shift', content):
                valid_cs_file = True
                file_type = "NMR-STAR"
            if '_Gen_dist_constraint' in lp_list and \
                    not self.is_loop_empty(input_data, '_Gen_dist_constraint', content):
                valid_re_file = True
                file_type = "NMR-STAR"
        else:
            err_msg = "Problem with input file"
            print(err_msg)
        return file_type, valid_cs_file, valid_re_file

    @staticmethod
    def is_loop_empty(stardata, lpcategory, content):
        is_empty = False
        if content == "Entry":
            lp_data = stardata.get_loops_by_category(lpcategory)
            for lpd in lp_data:
                if len(lpd.data) == 0:
                    is_empty = True
        elif content == "Saveframe":
            lp_data = stardata.get_loops_by_category(lpcategory)
            if len(lp_data.data) == 0:
                is_empty = True
        else:
            if len(stardata.data) == 0:
                is_empty = True
        return is_empty

    @staticmethod
    def get_sequence_from_nef(nefdata,
                              lp_category='nef_chemical_shift',
                              seq_id='sequence_code',
                              res_id='residue_name',
                              chain_id='chain_code'):
        """
        Redurns sequence from chemical shift loop for a given STAR data object
        :param nefdata: STAR data object
        :param lp_category: loop category
        :param seq_id: sequence number tag
        :param res_id: residue name tag
        :param chain_id: chain identifier tag
        :return: sequence
        """
        try:
            loop_data = nefdata.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loop_data = nefdata.get_loop_by_category(lp_category)
            except AttributeError:
                loop_data = [nefdata]
        seq = []

        for lp in loop_data:
            seq_dict = {}
            seqdat = lp.get_data_by_tag([seq_id, res_id, chain_id])
            chains = set([i[2] for i in seqdat])
            seq1 = sorted(set(['{}-{:03d}-{}'.format(i[2], int(i[0]), i[1]) for i in seqdat]))
            if len(seq1[0].split("-")[-1]) > 1:
                if len(chains) > 1:
                    for c in chains:
                        # seq2 = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq1 if i.split("-")[0] == c])
                        seq2 = [i.split("-")[-1] for i in seq1 if i.split("-")[0] == c]
                        seq_dict[c] = seq2
                else:
                    # seq2 = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq1])
                    seq2 = [i.split("-")[-1] for i in seq1]
                    seq_dict[list(chains)[0]] = seq2
            else:
                if len(chains) > 1:
                    for c in chains:
                        # seq2 = "".join([i.split("-")[-1] for i in seq1 if i.split("-")[0] == c])
                        seq2 = [i.split("-")[-1] for i in seq1 if i.split("-")[0] == c]
                        seq_dict[c] = seq2
                else:
                    # seq2 = "".join([i.split("-")[-1] for i in seq1])
                    seq2 = [i.split("-")[-1] for i in seq1]
                    seq_dict[list(chains)[0]] = seq2
            seq.append(seq_dict)
        return seq

    @staticmethod
    def get_sequence_from_nmrstar(stardata,
                                  lp_category='Atom_chem_shift',
                                  seq_id='Comp_index_ID',
                                  res_id='Comp_ID',
                                  chain_id='Entity_assembly_id'):
        """
        Returns sequence from chemical shift loop
        :param stardata: STAR data object
        :param lp_category: loop category
        :param seq_id: sequence ID tag
        :param res_id: residue name tag
        :param chain_id: chain ID tag
        :return: sequence
        """
        try:
            loop_data = stardata.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loop_data = stardata.get_loop_by_category(lp_category)
            except AttributeError:
                loop_data = [stardata]
        seq = []
        for csl in loop_data:
            seq_dict = {}
            if '_Atom_chem_shift.Entity_assembly_ID' not in csl.get_tag_names():
                seqdat = csl.get_data_by_tag([seq_id, res_id])
                for i in seqdat:
                    i.append(".")
            else:
                seqdat = csl.get_data_by_tag([seq_id, res_id, chain_id])

            chains = set([i[2] for i in seqdat])
            seq1 = sorted(set(['{}-{:03d}-{}'.format(i[2], int(i[0]), i[1]) for i in seqdat]))
            if len(seq1[0].split("-")[-1]) > 1:
                if len(chains) > 1:
                    for c in chains:
                        # seq2 = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq1 if i.split("-")[0] == c])
                        seq2 = [i.split("-")[-1] for i in seq1 if i.split("-")[0] == c]
                        seq_dict[c] = seq2
                else:
                    # seq2 = "".join([self.getOneLetter(i.split("-")[-1]) for i in seq1])
                    seq2 = [i.split("-")[-1] for i in seq1]
                    seq_dict[list(chains)[0]] = seq2
            else:
                if len(chains) > 1:
                    for c in chains:
                        # seq2 = "".join([i.split("-")[-1] for i in seq1 if i.split("-")[0] == c])
                        seq2 = [i.split("-")[-1] for i in seq1 if i.split("-")[0] == c]
                        seq_dict[c] = seq2
                else:
                    # seq2 = "".join([i.split("-")[-1] for i in seq1])
                    seq2 = [i.split("-")[-1] for i in seq1]
                    seq_dict[list(chains)[0]] = seq2
            seq.append(seq_dict)
        return seq

    @staticmethod
    def get_saveframes_and_loops(stardata, content_type):
        """
        Returns the list of saveframe categories and loop categories for a given stardata object
        :param stardata: stardata object
        :param content_type: Entry/Saveframe/Loop
        :return: list of saveframe categories, loop categories
        """
        sf_list = []
        lp_list = []
        if content_type == "Entry":
            for sf in stardata.frame_list:
                sf_list.append(sf.category)
                for lp in sf:
                    lp_list.append(lp.category)
        elif content_type == "Saveframe":
            for lp in stardata:
                lp_list.append(lp.category)
        else:
            lp_list.append(stardata.category)
        return sf_list, lp_list

    def validate_atom_nomenclature(self, stardata,
                                   lp_category='Atom_chem_shift',
                                   seq_id='Comp_index_ID',
                                   res_id='Comp_ID',
                                   atom_id='Atom_ID'):
        """
        Validates the atoms in a given loop against IUPAC standard
        :param stardata: star data object
        :param lp_category: loop category
        :param seq_id: sequence id tag
        :param res_id: residue id tag
        :param atom_id: atom id tag
        :return: list of non standard atoms
        """

        try:
            loop_data = stardata.get_loops_by_category(lp_category)
        except AttributeError:
            try:
                loop_data = [stardata.get_loop_by_category(lp_category)]
            except AttributeError:
                loop_data = [stardata]

        ns = []
        for lp in loop_data:
            try:
                atm_data = lp.get_data_by_tag([seq_id, res_id, atom_id])
                for i in atm_data:
                    try:
                        if i[2] not in self.atom_dict[i[1].upper()]:
                            ns.append(i)
                    except KeyError:
                        ns.append(i)
            except ValueError:
                print("One of the following tag is missing ", seq_id, res_id, atom_id)

            # nonStandard = [i for i in atm_data if i[2] not in self.atomDict[i[1].upper]]
            # ns.append(nonStandard)
        return ns

    def get_nmrstar_tag(self, tag):
        n = self.tag_map[0].index(tag)
        return [self.tag_map[1][n], self.tag_map[2][n]]

    def get_nef_tag(self, tag):
        n = self.tag_map[1].index(tag)
        return self.tag_map[0][n]

    def get_nmrstar_atom(self, res, nef_atom):
        atom = None
        ambiguity_code = 1
        try:
            atms = self.atom_dict[res]
            atom_list = []
            try:
                refatm = re.findall(r'(\S+)([xyXY])([%*])$|(\S+)([%*])$|(\S+)([xyXY]$)', nef_atom)[0]
                atm_set = [refatm.index(i) for i in refatm if i != ""]
                if atm_set == [0, 1, 2]:
                    atom = refatm[0]
                    pattern = re.compile(r'%s\S\d+' % (refatm[0]))
                    alist2 = [i for i in atms if re.search(pattern, i)]
                    xid = sorted(set([int(i[len(refatm[0])]) for i in alist2]))
                    if refatm[1] == "x" or refatm[1] == "X":
                        atom_list = [i for i in alist2 if int(i[len(refatm[0])]) == xid[0]]
                    else:
                        atom_list = [i for i in alist2 if int(i[len(refatm[0])]) == xid[1]]
                    ambiguity_code = 2
                elif atm_set == [3, 4]:
                    atom = refatm[3]
                    if refatm[4] == "%":
                        pattern = re.compile(r'%s\d+' % (refatm[3]))
                    elif refatm[4] == "*":
                        pattern = re.compile(r'%s\S+' % (refatm[3]))
                    else:
                        pattern = None
                        print("Wrong NEF atom {}".format(nef_atom))
                    atom_list = [i for i in atms if re.search(pattern, i)]
                    ambiguity_code = 1

                elif atm_set == [5, 6]:
                    atom = refatm[5]
                    pattern = re.compile(r'%s\S+' % (refatm[5]))
                    atom_list = [i for i in atms if re.search(pattern, i)]
                    if len(atom_list) != 2:
                        atom_list = []
                    elif refatm[6] == "y" or refatm[6] == "Y":
                        # alist.reverse()[]
                        atom_list = atom_list[-1:]
                    elif refatm[6] == "x" or refatm[6] == "X":
                        atom_list = atom_list[:1]
                    else:
                        print("Wrong NEF atom {}".format(nef_atom))
                    ambiguity_code = 2

                else:
                    print("Wrong NEF atom {}".format(nef_atom))
            except IndexError:

                # print nefAtom
                pass
                atom = nef_atom
            if len(atom_list) == 0:
                if nef_atom in atms:
                    atom_list.append(nef_atom)
                else:
                    if nef_atom == "H%":  # To handle terminal protons
                        atom_list = ['H1', 'H2', 'H3']
                        atom = "H"
        except KeyError:
            # self.logfile.write("%s\tResidue not found,%s,%s\n"%(self.TimeStamp(time.time()),res,nefAtom))
            # print "Residue not found",res,nefAtom
            if res != ".":
                print("Non-standard residue found {}".format(res))
            atom_list = []
            atom = nef_atom

            if nef_atom == "H%":
                atom_list = ['H1', 'H2', 'H3']
                atom = "H"
        return [atom, atom_list, ambiguity_code]

    def get_nmrstar_loop_tags(self, neflooptags):
        aut_tag = []
        nt = []
        for t in neflooptags:
            st = self.get_nmrstar_tag(t)
            if st[0] != st[1]:
                aut_tag.append(st[1])
            nt.append(st[0])
        if len(aut_tag) != 0:
            out_tag = nt + aut_tag
        else:
            out_tag = nt
        if neflooptags[0].split(".")[0] == "_nef_chemical_shift":
            out_tag.append('_Atom_chem_shift.Ambiguity_code')
            out_tag.append('_Atom_chem_shift.Ambiguity_set_ID')
            out_tag.append('_Atom_chem_shift.Assigned_chem_shift_list_ID')
        if neflooptags[0].split(".")[0] == "_nef_distance_restraint":
            out_tag.append('_Gen_dist_constraint.Member_logic_code')
        return out_tag

    @staticmethod
    def time_stamp(ts):
        return datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

    def translate_row(self, f_tags, t_tags, row_data):
        # print (f_tags)
        out_row = []
        res_list = self.get_identifier(f_tags)
        # print (res_list)
        tmp_dict = {}
        for res1 in res_list:
            try:
                tmp_dict[res1[0]] = self.seq_dict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][0]
            except KeyError:
                tmp_dict[res1[0]] = row_data[f_tags.index(res1[0])]
            try:
                tmp_dict[res1[1]] = self.seq_dict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][1]
            except KeyError:
                tmp_dict[res1[1]] = row_data[f_tags.index(res1[1])]
        # print (tmp_dict)
        if len(f_tags) != len(t_tags):
            out = [None] * len(t_tags)
            for j in f_tags:
                stgs = self.get_nmrstar_tag(j)
                if stgs[0] == stgs[1]:
                    out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                else:
                    if 'chain_code' in j or 'sequence_code' in j:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = tmp_dict[j]
                    else:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]

                # else:
                #   print ("ERROR",f_tags)
            out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def translate_seq_row(self, f_tags, t_tags, row_data):
        out_row = []
        if len(f_tags) != len(t_tags):
            out = [None] * len(t_tags)
            for j in f_tags:
                stgs = self.get_nmrstar_tag(j)
                if stgs[0] == stgs[1]:
                    out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                else:
                    if j == '_nef_sequence.chain_code':
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = self.chains.index(row_data[f_tags.index(j)]) + 1
                    elif j == '_nef_sequence.sequence_code':
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = \
                            self.cid[self.chains.index(row_data[f_tags.index('_nef_sequence.chain_code')])]
                    else:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]
            out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def translate_cs_row(self, f_tags, t_tags, row_data):
        new_id = None
        out_row = []

        if '_nef_chemical_shift.chain_code' in f_tags and '_nef_chemical_shift.sequence_code' in f_tags:

            cci = f_tags.index('_nef_chemical_shift.chain_code')
            sci = f_tags.index('_nef_chemical_shift.sequence_code')
            try:
                old_id = [i for i in self.seq_dict.keys() if i[0] == row_data[cci] and i[1] == row_data[sci]][0]
                new_id = self.seq_dict[old_id]
            except AttributeError:
                new_id = (cci, sci)
            except IndexError:
                new_id = (cci, sci)
        if len(f_tags) != len(t_tags):
            atm_index = f_tags.index('_nef_chemical_shift.atom_name')
            res_index = f_tags.index('_nef_chemical_shift.residue_name')
            n_atm = self.get_nmrstar_atom(row_data[res_index], row_data[atm_index])[1]
            ambi = self.get_nmrstar_atom(row_data[res_index], row_data[atm_index])[2]

            for i in n_atm:
                out = [None] * len(t_tags)
                for j in f_tags:
                    stgs = self.get_nmrstar_tag(j)
                    if stgs[0] == stgs[1]:
                        out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                    else:
                        if j == '_nef_chemical_shift.atom_name':
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = i
                        elif j == '_nef_chemical_shift.chain_code':
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = new_id[0]
                        elif j == '_nef_chemical_shift.sequence_code':
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = new_id[1]
                        else:
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                            out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]
                    out[t_tags.index('_Atom_chem_shift.Ambiguity_code')] = ambi
                    out[t_tags.index('_Atom_chem_shift.Ambiguity_set_ID')] = '.'
                out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    @staticmethod
    def get_identifier(tag_list):
        out_list = []  # type: List[]
        for j in range(1, 16):
            out = [None] * 2
            chk_string = re.compile('\S+.chain_code_{}'.format(j))
            r1 = [chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
            if len(r1) > 0:
                out[0] = r1[0]
            chk_string = re.compile('\S+.sequence_code_{}'.format(j))
            r2 = [chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
            if len(r2) > 0:
                out[1] = r2[0]
            if len(r1) > 0 and len(r2) > 0:
                out_list.append(out)
        #             chk_string = re.compile('\S+.residue_name_{}'.format(j))
        #             r=[chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
        #             if len(r)>0: out[2]=r[0]
        #             chk_string = re.compile('\S+.atom_name_{}'.format(j))
        #             r=[chk_string.search(i).group() for i in tag_list if chk_string.search(i)]
        #             if len(r)>0: out[3]=r[0]
        return out_list

    def translate_restraint_row(self, f_tags, t_tags, row_data):
        out_row = []
        res_list = self.get_identifier(f_tags)
        # print (res_list)
        tmp_dict = {}
        for res1 in res_list:
            try:
                tmp_dict[res1[0]] = self.seq_dict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][0]
            except KeyError:
                tmp_dict[res1[0]] = row_data[f_tags.index(res1[0])]
            try:
                tmp_dict[res1[1]] = self.seq_dict[(row_data[f_tags.index(res1[0])], row_data[f_tags.index(res1[1])])][1]
            except KeyError:
                tmp_dict[res1[1]] = row_data[f_tags.index(res1[1])]
        if len(f_tags) != len(t_tags):
            atm_index1 = f_tags.index('_nef_distance_restraint.atom_name_1')
            res_index1 = f_tags.index('_nef_distance_restraint.residue_name_1')
            atm_index2 = f_tags.index('_nef_distance_restraint.atom_name_2')
            res_index2 = f_tags.index('_nef_distance_restraint.residue_name_2')
            n_atm1 = self.get_nmrstar_atom(row_data[res_index1], row_data[atm_index1])[1]
            n_atm2 = self.get_nmrstar_atom(row_data[res_index2], row_data[atm_index2])[1]

            for i in n_atm1:
                for k in n_atm2:
                    out = [None] * len(t_tags)
                    for j in f_tags:
                        stgs = self.get_nmrstar_tag(j)
                        if stgs[0] == stgs[1]:
                            out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                        else:
                            if j == '_nef_distance_restraint.atom_name_1':
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = i
                            elif 'chain_code_1' in j or 'sequence_code_1' in j:
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = tmp_dict[j]
                            elif j == '_nef_distance_restraint.atom_name_2':
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = k
                            elif 'chain_code_2' in j or 'sequence_code_2' in j:
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = tmp_dict[j]
                            else:
                                out[t_tags.index(stgs[0])] = row_data[f_tags.index(j)]
                                out[t_tags.index(stgs[1])] = row_data[f_tags.index(j)]

                    out_row.append(out)
        else:
            out_row.append(row_data)
        return out_row

    def nef_to_nmrstar(self, nef_file):
        (filePath, fileName) = ntpath.split(os.path.realpath(nef_file))
        is_done = True
        info = []
        warning = []
        error = []
        star_file = filePath + "/" + fileName.split(".")[0] + ".str"
        (isReadable, dat_content, nefData) = self.read_input(nef_file)
        try:
            star_data = pynmrstar.Entry.from_scratch(nefData.entry_id)
        except AttributeError:
            star_data = pynmrstar.Entry.from_scratch(fileName.split(".")[0])
            warning.append('Not a complete Entry')
        if isReadable:
            if dat_content == "Entry":
                self.chains = sorted(list(set(nefData.get_loops_by_category('nef_sequence')[0].get_tag('chain_code'))))
            elif dat_content == "Saveframe":
                self.chains = sorted(list(set(nefData[0].get_tag('chain_code'))))
            elif dat_content == "Loop":
                self.chains = sorted(list(set(nefData.get_tag('chain_code'))))
            else:
                is_done = False
                error.append('File content unknown')

            cs_list = 0
            if dat_content == "Entry":
                for saveframe in nefData:
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)

                    for tag in saveframe.tags:
                        if tag[0].lower() == "sf_category":
                            sf.add_tag("Sf_category", self.get_nmrstar_tag(saveframe.category)[0])
                        else:
                            neftag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            sf.add_tag(self.get_nmrstar_tag(neftag)[0], tag[1])
                    if saveframe.category == "nef_nmr_meta_data":
                        sf.add_tag("NMR_STAR_version", "3.2.0.15")
                        sf.add_tag("Generated_date", self.time_stamp(time.time()), update=True)
                    for loop in saveframe:
                        if loop.category == "_nef_sequence":
                            self.cid = []  # Comp_index_ID list
                            for c in self.chains:  # Comp_index_ID initialized with 1
                                self.cid.append(1)
                            self.seq_dict = {}

                        if loop.category == '_nef_distance_restraint':
                            r_index_id = 1
                        if loop.category == "_nef_chemical_shift":
                            cs_list += 1
                        lp = pynmrstar.Loop.from_scratch()
                        lp_cols = self.get_nmrstar_loop_tags(loop.get_tag_names())
                        for t in lp_cols:
                            lp.add_tag(t)
                        # print (loop.category,lp.category,lp.get_tag_names(),loop.get_tag_names())
                        for dat in loop.data:
                            if loop.category == "_nef_sequence":
                                dd = self.translate_seq_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                self.cid[
                                    self.chains.index(dat[loop.get_tag_names().index('_nef_sequence.chain_code')])] += 1
                                for d in dd:
                                    lp.add_data(d)
                                    self.seq_dict[(dat[loop.get_tag_names().index('_nef_sequence.chain_code')],
                                                   dat[loop.get_tag_names().index('_nef_sequence.sequence_code')])] = (
                                        d[lp.get_tag_names().index('_Chem_comp_assembly.Entity_assembly_ID')],
                                        d[lp.get_tag_names().index('_Chem_comp_assembly.Comp_index_ID')])

                            elif loop.category == "_nef_chemical_shift":
                                dd = self.translate_cs_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    d[lp.get_tag_names().index(
                                        '_Atom_chem_shift.Assigned_chem_shift_list_ID')] = cs_list
                                    lp.add_data(d)
                            elif loop.category == '_nef_distance_restraint':
                                dd = self.translate_restraint_row(loop.get_tag_names(), lp.get_tag_names(), dat)

                                for d in dd:
                                    d[lp.get_tag_names().index('_Gen_dist_constraint.Index_ID')] = r_index_id
                                    if len(dd) > 1:
                                        d[lp.get_tag_names().index('_Gen_dist_constraint.Member_logic_code')] = "OR"
                                    lp.add_data(d)
                                    r_index_id += 1
                            else:
                                dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                                for d in dd:
                                    lp.add_data(d)

                        # print (loop.data[0])
                        sf.add_loop(lp)
                    star_data.add_saveframe(sf)
                star_data.normalize()
                with open(star_file, 'w') as wstarfile:
                    wstarfile.write(str(star_data))
            elif dat_content == "Saveframe" or dat_content == "Loop":
                if dat_content == "Saveframe":
                    saveframe = nefData
                    sf = pynmrstar.Saveframe.from_scratch(saveframe.name)
                    for tag in saveframe.tags:
                        if tag[0].lower() == "sf_category":
                            try:
                                sf.add_tag("Sf_category", self.get_nmrstar_tag(saveframe.category)[0])
                            except ValueError:
                                sf.add_tag("Sf_category", self.get_nmrstar_tag(tag[1])[0])
                        else:
                            neftag = '{}.{}'.format(saveframe.tag_prefix, tag[0])
                            sf.add_tag(self.get_nmrstar_tag(neftag)[0], tag[1])
                    if saveframe.category == "nef_nmr_meta_data":
                        sf.add_tag("NMR_STAR_version", "3.2.0.15")

                else:
                    sf = pynmrstar.Saveframe.from_scratch(nefData.category)
                    if nefData.category == "_nef_chemical_shift":
                        sf.add_tag("_Assigned_chem_shift_list.Sf_category", 'nef_chemical_shift')
                    saveframe = [nefData]
                for loop in saveframe:
                    if loop.category == "_nef_sequence":
                        self.cid = []  # Comp_index_ID list
                        for c in self.chains:  # Comp_index_ID initialized with 1
                            self.cid.append(1)
                        self.seq_dict = {}

                    if loop.category == '_nef_distance_restraint':
                        r_index_id = 1
                    if loop.category == "_nef_chemical_shift":
                        cs_list += 1
                    lp = pynmrstar.Loop.from_scratch()
                    lp_cols = self.get_nmrstar_loop_tags(loop.get_tag_names())
                    for t in lp_cols:
                        lp.add_tag(t)
                    # print (loop.category,lp.category,lp.get_tag_names(),loop.get_tag_names())
                    for dat in loop.data:
                        if loop.category == "_nef_sequence":
                            dd = self.translate_seq_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            self.cid[
                                self.chains.index(dat[loop.get_tag_names().index('_nef_sequence.chain_code')])] += 1
                            for d in dd:
                                lp.add_data(d)
                                self.seq_dict[(dat[loop.get_tag_names().index('_nef_sequence.chain_code')],
                                               dat[loop.get_tag_names().index('_nef_sequence.sequence_code')])] = (
                                    d[lp.get_tag_names().index('_Chem_comp_assembly.Entity_assembly_ID')],
                                    d[lp.get_tag_names().index('_Chem_comp_assembly.Comp_index_ID')])

                        elif loop.category == "_nef_chemical_shift":
                            dd = self.translate_cs_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            for d in dd:
                                d[lp.get_tag_names().index('_Atom_chem_shift.Assigned_chem_shift_list_ID')] = cs_list
                                lp.add_data(d)
                        elif loop.category == '_nef_distance_restraint':
                            dd = self.translate_restraint_row(loop.get_tag_names(), lp.get_tag_names(), dat)

                            for d in dd:
                                d[lp.get_tag_names().index('_Gen_dist_constraint.Index_ID')] = r_index_id
                                if len(dd) > 1:
                                    d[lp.get_tag_names().index('_Gen_dist_constraint.Member_logic_code')] = "OR"
                                lp.add_data(d)
                                r_index_id += 1
                        else:
                            dd = self.translate_row(loop.get_tag_names(), lp.get_tag_names(), dat)
                            for d in dd:
                                lp.add_data(d)

                    # print (loop.data[0])
                    sf.add_loop(lp)
                star_data.add_saveframe(sf)
            star_data.normalize()
            with open(star_file, 'w') as wstarfile:
                wstarfile.write(str(star_data))

        else:
            is_done = False
            error.append('Input file not readable')
        return [is_done, json.dumps({'INFO': info, 'WARNING': warning, 'ERROR': error})]
