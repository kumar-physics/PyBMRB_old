#!/usr/bin/env python

from unittest import TestCase

from Translator import Translator


class TestTranslator(TestCase):
    def setUp(self):
        self.Translator = Translator()

    def test_read_input(self):
        file_readble, content_type, data = self.Translator.read_input("test/1nk2.nef")
        self.assertTrue(file_readble)
        self.assertEqual(content_type, "Entry")
        file_readble, content_type, data = self.Translator.read_input("test/1nk2.str")
        self.assertTrue(file_readble)
        self.assertEqual(content_type, "Entry")
        file_readble, content_type, data = self.Translator.read_input("test/unit_cs_sf.nef")
        self.assertTrue(file_readble)
        self.assertEqual(content_type, "Saveframe")
        file_readble, content_type, data = self.Translator.read_input("test/unit_cs_sf.str")
        self.assertTrue(file_readble)
        self.assertEqual(content_type, "Saveframe")
        file_readble, content_type, data = self.Translator.read_input("test/unit_cs_sf2.str")
        self.assertFalse(file_readble)
        self.assertEqual(content_type, None)

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

    def test_load_map_file(self):
        self.Translator.load_map_file()
        self.assertEqual([len(i) for i in self.Translator.tag_map], [258, 258, 258],
                         msg="There are 258 tags defined in the lib file")

    def test_load_nef_info(self):
        self.Translator.load_nef_info()
        self.assertEqual(len(self.Translator.nef_info), 270,
                         msg="There are 270 mandatory tags for NEF defined in the file")

    # def test_validate_file(self):
    #     self.fail()
    #
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
        self.assertEqual(seq_data[0].keys(), ['A'])
        self.assertEqual(len(seq_data[0]['A']), 214)
        self.assertEqual(seq_data[1].keys(), ['B'])
        self.assertEqual(len(seq_data[1]['B']), 5)
        self.assertEqual(seq_data[2].keys(), ['C'])
        self.assertEqual(len(seq_data[2]['C']), 5)

    def test_get_sequence_from_nmrstar(self):
        readable, content_type, stardata = self.Translator.read_input('test/test_seq.str')
        seq_data = self.Translator.get_sequence_from_nmrstar(stardata)
        self.assertEqual(len(seq_data), 3)
        self.assertEqual(seq_data[0].keys(), ['1'])
        self.assertEqual(len(seq_data[0]['1']), 214)
        self.assertEqual(seq_data[1].keys(), ['2'])
        self.assertEqual(len(seq_data[1]['2']), 5)
        self.assertEqual(seq_data[2].keys(), ['3'])
        self.assertEqual(len(seq_data[2]['3']), 5)

    # def test_get_saveframes_and_loops(self):
    #     self.fail()
    #
    # def test_validate_atom_nomenclature(self):
    #     self.fail()
    #
    # def test_get_nmrstar_tag(self):
    #     self.fail()
    #
    # def test_get_nef_tag(self):
    #     self.fail()
    #
    # def test_get_nmrstar_atom(self):
    #     self.fail()
    #
    # def test_get_nmrstar_loop_tags(self):
    #     self.fail()
    #
    # def test_time_stamp(self):
    #     self.fail()
    #
    # def test_translate_row(self):
    #     self.fail()
    #
    # def test_translate_seq_row(self):
    #     self.fail()
    #
    # def test_translate_cs_row(self):
    #     self.fail()
    #
    # def test_get_identifier(self):
    #     self.fail()
    #
    # def test_translate_restraint_row(self):
    #     self.fail()
    #
    # def test_nef_to_nmrstar(self):
    #     self.fail()
