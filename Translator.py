from __future__ import print_function

import csv
import json
import sys
import os
import ntpath

PY3 = (sys.version_info[0] == 3)

(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))

try:
    import pynmrstar

except ImportError as e:
    sys.path.append(scriptPath + '/PyNMRSTAR')
    try:
        import bmrb as pynmrstar
    except ImportError as e:
        print(scriptPath, scriptName)
        print("ERROR: PyNMRSTAR parser is not available")
        print(str(e))
        exit(1)

__version__ = "v1.0"


class Translator:
    map_file = scriptPath + '/lib/NEF_NMRSTAR_equivalence.csv'
    nef_info_file = scriptPath + '/lib/NEF_mandatory,csv'
    atom_nomenclature_file = scriptPath + '/lib/atomDict.json'
    amino_acids_code_file = scriptPath + '/lib/codeDict.json'

    def __init__(self):
        self.atom_dict = None
        self.code_dict = None
        self.tag_map = None
        self.nef_info = None
        self.load_atom_dict()
        self.load_code_dict()

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
        """Reads the codeDict.json file and creates a dictionary of residues and atoms"""
        try:
            with open(self.amino_acids_code_file, 'r') as code_file:
                self.code_dict = json.loads(code_file.read())
        except IOError:
            err_msg = "codeDict.json file is missing! Check whether the file is inside {}".format(scriptPath + '/lib')
            print(err_msg)

    def get_one_letter_code(self, res):
        """
        Takes 3 letter code as input and return one letter code. If one letter code not defined in the dictionary
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
        try:
            with open(self.map_file, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                map_data = []
                for r in spamreader:
                    if r[0][0] != '#':
                        map_data.append(r)
            self.tag_map = list(map(list, zip(*map_data)))
        except IOError:
            err_msg = "NEF-NMRSTAR_equivalence.csv file is missing! check the file is inside  {} ".format(
                scriptPath + '/lib')
            print(err_msg)
            self.tag_map = []

    def load_nef_info(self):
        """
        Reads mandatory tag information for NEF file
        """
        try:
            with open(self.nef_info_file, 'r') as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                map_data = []
                for r in spamreader:
                    if r[0][0] != '#':
                        map_data.append(r)
            self.nef_info = map_data
        except IOError:
            err_msg = "NEF_mandatory.csv file is missing! check the file is inside  {} ".format(scriptPath + '/lib')
            print(err_msg)
            self.nef_info = []

    def validate_file(self, input_file):
        valid_file = False
        readable, content, input_data = self.read_input(input_file)
        if readable:
            pass
        else:
            err_msg = "Problem with input file"
            print(err_msg)

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
