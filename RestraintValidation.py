#!/usr/bin/env python

from __future__ import print_function

import sys
import os
import ntpath


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
        self.cfile = 'test/data/1nk2.cif'
        self.rfile = None
        self.cif_data = None
        self.max_models = None
        self.read_ciffile()
        pass

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


if __name__ == "__main__":
    p = RestraintValidation()
