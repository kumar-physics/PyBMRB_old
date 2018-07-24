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
        self.assertEqual([len(i) for i in self.Translator.tag_map], [258,258,258],
                         msg="There are 258 tags defined in the lib file")

    def test_load_nef_info(self):
        self.Translator.load_nef_info()
        self.assertEqual(len(self.Translator.nef_info), 270,
                         msg="There are 270 mandatory tags for NEF defined in the file")

    # def test_validate_file(self):
    #     self.fail()
    #
    # def test_is_loop_empty(self):
    #     self.fail()
    #
    # def test_get_sequence_from_nef(self):
    #     self.fail()
    #
    # def test_get_sequence_from_nmrstar(self):
    #     self.fail()
    #
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
