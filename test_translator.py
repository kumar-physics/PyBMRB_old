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

    # def test_load_atom_dict(self):
    #     self.fail()
    #
    # def test_load_code_dict(self):
    #     self.fail()
    #
    # def test_get_one_letter_code(self):
    #     self.fail()
    #
    # def test_load_map_file(self):
    #     self.fail()
    #
    # def test_load_nef_info(self):
    #     self.fail()
    #
    # def test_validate_file(self):
    #     self.fail()
    #
    # def test_get_saveframes_and_loops(self):
    #     self.fail()
