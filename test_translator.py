#!/usr/bin/env python

from unittest import TestCase

from Translator import Translator


class TestTranslator(TestCase):
    def setUp(self):
        self.Translator = Translator()

    def test_read_input(self):
        file_readble, content_type, data = self.Translator.read_input("test/test_good_entry.nef")
        self.assertTrue(file_readble)
        self.assertAlmostEqual(content_type, "Entry")
        file_readble, content_type, data = self.Translator.read_input("test/test_good_saveframe.nef")
        self.assertTrue(file_readble)
        self.assertAlmostEqual(content_type, "Saveframe")
        file_readble, content_type, data = self.Translator.read_input("test/test_good_loop.nef")
        self.assertTrue(file_readble)
        self.assertAlmostEqual(content_type, "Loop")
        file_readble, content_type, data = self.Translator.read_input("test/test_good_entry.str")
        self.assertTrue(file_readble)
        self.assertAlmostEqual(content_type, "Entry")
        file_readble, content_type, data = self.Translator.read_input("test/test_good_saveframe.str")
        self.assertTrue(file_readble)
        self.assertAlmostEqual(content_type, "Saveframe")
        file_readble, content_type, data = self.Translator.read_input("test/test_good_loop.str")
        self.assertTrue(file_readble)
        self.assertAlmostEqual(content_type, "Loop")
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_entry.nef")
        self.assertFalse(file_readble)
        self.assertIsNone(content_type)
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_saveframe.nef")
        self.assertFalse(file_readble)
        self.assertIsNone(content_type)
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_loop.nef")
        self.assertFalse(file_readble)
        self.assertIsNone(content_type)
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_entry.str")
        self.assertFalse(file_readble)
        self.assertIsNone(content_type)
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_saveframe.str")
        self.assertFalse(file_readble)
        self.assertIsNone(content_type)
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_loop.str")
        self.assertFalse(file_readble)
        self.assertIsNone(content_type)

    def test_load_atom_dict(self):
        self.Translator.load_atom_dict()
        self.assertEqual(len(self.Translator.atom_dict), 37,
                         msg="There should be 37 standard residues defined in the dictionary")

    def test_load_code_dict(self):
        self.Translator.load_code_dict()
        self.assertEqual(len(self.Translator.code_dict), 37,
                         msg="There should be 37 standard residues defined in the dictionary")

    def test_get_one_letter_code(self):
        self.assertEqual(self.Translator.get_one_letter_code('ALA'), 'A')
        self.assertEqual(self.Translator.get_one_letter_code('val'), 'V')
        self.assertEqual(self.Translator.get_one_letter_code('mno'), '?')
        self.assertEqual(self.Translator.get_one_letter_code('DA'), 'A')

    def test_load_map_file(self):
        self.Translator.load_map_file()
        self.assertEqual([len(i) for i in self.Translator.tag_map], [258, 258, 258],
                         msg="There are 258 tags defined in the lib file")

    def test_load_nef_info(self):
        self.Translator.load_nef_info()
        self.assertEqual(len(self.Translator.nef_info), 270,
                         msg="There are 270 mandatory tags for NEF defined in the file")

    def test_validate_file(self):
        file_type, cs_info, rt_info = self.Translator.validate_file('test/test_good_entry.nef')
        self.assertEqual(file_type, 'NEF')
        self.assertTrue(cs_info)
        self.assertTrue(rt_info)
        file_type, cs_info, rt_info = self.Translator.validate_file('test/test_good_entry.str')
        self.assertEqual(file_type, 'NMR-STAR')
        self.assertTrue(cs_info)
        self.assertTrue(rt_info)
        file_type, cs_info, rt_info = self.Translator.validate_file('test/test_good_loop.str')
        self.assertEqual(file_type, 'NMR-STAR')
        self.assertTrue(cs_info)
        self.assertFalse(rt_info)
        file_type, cs_info, rt_info = self.Translator.validate_file('test/test_good_loop.nef')
        self.assertEqual(file_type, 'NEF')
        self.assertTrue(cs_info)
        self.assertFalse(rt_info)
        file_type, cs_info, rt_info = self.Translator.validate_file('test/test_bad_entry.str')
        self.assertIsNone(file_type)
        self.assertFalse(cs_info)
        self.assertFalse(rt_info)
        file_type, cs_info, rt_info = self.Translator.validate_file('test/test_bad_saveframe.nef')
        self.assertIsNone(file_type)
        self.assertFalse(cs_info)
        self.assertFalse(rt_info)
        file_type, cs_info, rt_info = self.Translator.validate_file('test/loop_test.nef')
        self.assertEqual(file_type, 'NEF')
        self.assertFalse(cs_info)
        self.assertTrue(rt_info)

    def test_is_loop_empty(self):
        readable, content_type, stardata = self.Translator.read_input('test/loop_test.nef')
        self.assertFalse(self.Translator.is_loop_empty(stardata, '_nef_sequence', content_type))
        self.assertTrue(self.Translator.is_loop_empty(stardata, '_nef_chemical_shift', content_type))
        readable, content_type, stardata = self.Translator.read_input('test/loop_test.str')
        self.assertTrue(self.Translator.is_loop_empty(stardata, '_Chem_comp_assembly', content_type))
        self.assertFalse(self.Translator.is_loop_empty(stardata, '_Atom_chem_shift', content_type))

    def test_get_sequence_from_nef(self):
        readable, content_type, nefdata = self.Translator.read_input('test/test_seq.nef')
        seq_data = self.Translator.get_sequence_from_nef(nefdata)
        self.assertEqual(len(seq_data), 3)
        self.assertEqual(list(seq_data[0].keys()), ['A'])
        self.assertEqual(len(seq_data[0]['A']), 214)
        self.assertEqual(list(seq_data[1].keys()), ['B'])
        self.assertEqual(len(seq_data[1]['B']), 5)
        self.assertEqual(list(seq_data[2].keys()), ['C'])
        self.assertEqual(len(seq_data[2]['C']), 5)

    def test_get_sequence_from_nmrstar(self):
        readable, content_type, stardata = self.Translator.read_input('test/test_seq.str')
        seq_data = self.Translator.get_sequence_from_nmrstar(stardata)
        self.assertEqual(len(seq_data), 3)
        self.assertEqual(list(seq_data[0].keys()), ['1'])
        self.assertEqual(len(seq_data[0]['1']), 214)
        self.assertEqual(list(seq_data[1].keys()), ['2'])
        self.assertEqual(len(seq_data[1]['2']), 5)
        self.assertEqual(list(seq_data[2].keys()), ['3'])
        self.assertEqual(len(seq_data[2]['3']), 5)

    def test_get_saveframes_and_loops(self):
        readable, content_type, stardata = self.Translator.read_input('test/test_seq.str')
        sf_list, lp_list = self.Translator.get_saveframes_and_loops(stardata, content_type)
        self.assertEqual(sf_list, ['entry_information', 'assembly', 'assigned_chemical_shifts',
                                   'assigned_chemical_shifts', 'assigned_chemical_shifts',
                                   'general_distance_constraints', 'general_distance_constraints',
                                   'torsion_angle_constraints'])
        self.assertEqual(lp_list, ['_Software_applied_methods', '_Chem_comp_assembly', '_Atom_chem_shift',
                                   '_Atom_chem_shift', '_Atom_chem_shift', '_Gen_dist_constraint',
                                   '_Gen_dist_constraint', '_Torsion_angle_constraint'])
        readable, content_type, stardata = self.Translator.read_input('test/test_good_saveframe.str')
        sf_list, lp_list = self.Translator.get_saveframes_and_loops(stardata, content_type)
        self.assertEqual(sf_list, [])
        self.assertEqual(lp_list, ['_Atom_chem_shift'])

    def test_validate_atom_nomenclature(self):
        file_readble, content_type, data = self.Translator.read_input("test/test_bad_nomenclature.str")
        non_standard_list = self.Translator.validate_atom_nomenclature(data)
        self.assertEqual(len(non_standard_list), 17)
        file_readble, content_type, data = self.Translator.read_input("test/test_good_entry.str")
        non_standard_list = self.Translator.validate_atom_nomenclature(data)
        self.assertEqual(len(non_standard_list), 0)

    def test_get_nmrstar_tag(self):
        self.assertEqual(self.Translator.get_nmrstar_tag('_nef_sequence.sequence_code'),
                         ['_Chem_comp_assembly.Auth_seq_ID', '_Chem_comp_assembly.Comp_index_ID'])
        self.assertEqual(self.Translator.get_nmrstar_tag('_nef_chemical_shift.value'),
                         ['_Atom_chem_shift.Val', '_Atom_chem_shift.Val'])
        self.assertEqual(self.Translator.get_nmrstar_tag('nef_nmr_meta_data'),
                         ['entry_information', 'entry_information'])
        self.assertIsNone(self.Translator.get_nmrstar_tag('some_unknown_tag'))

    def test_get_nef_tag(self):
        self.assertEqual(self.Translator.get_nef_tag('_Chem_comp_assembly.Auth_seq_ID'),
                         '_nef_sequence.sequence_code')
        self.assertEqual(self.Translator.get_nef_tag('_Atom_chem_shift.Val'),
                         '_nef_chemical_shift.value')
        self.assertEqual(self.Translator.get_nef_tag('entry_information'),
                         'nef_nmr_meta_data')
        self.assertIsNone(self.Translator.get_nef_tag('some_unknown_tag'))

    def test_get_nmrstar_atom(self):
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('CYS', 'HB%')
        self.assertEqual(atom, 'HB')
        self.assertEqual(atom_list, ['HB2', 'HB3'])
        self.assertEqual(ambiguity_code,1)
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('TRP', 'CE%')
        self.assertEqual(atom, 'CE')
        self.assertEqual(atom_list, ['CE2', 'CE3'])
        self.assertEqual(ambiguity_code, 1)
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('TRP', 'CEX')
        self.assertEqual(atom, 'CE')
        self.assertEqual(atom_list, ['CE2'])
        self.assertEqual(ambiguity_code, 2)
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('TRP', 'CEy')
        self.assertEqual(atom, 'CE')
        self.assertEqual(atom_list, ['CE3'])
        self.assertEqual(ambiguity_code, 2)
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('LEU', 'HDY%')
        self.assertEqual(atom, 'HD')
        self.assertEqual(atom_list, ['HD21', 'HD22', 'HD23'])
        self.assertEqual(ambiguity_code, 2)
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('LEU', 'HD1%')
        self.assertEqual(atom, 'HD1')
        self.assertEqual(atom_list, ['HD11', 'HD12', 'HD13'])
        self.assertEqual(ambiguity_code, 1)
        atom, atom_list, ambiguity_code = self.Translator.get_nmrstar_atom('LEU', 'HD*')
        self.assertEqual(atom, 'HD')
        self.assertEqual(atom_list, ['HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'])
        self.assertEqual(ambiguity_code, 1)



    def test_nef_to_nmrstar(self):
        status, jdump = self.Translator.nef_to_nmrstar('test/test_translation.nef')
        self.assertTrue(status)
