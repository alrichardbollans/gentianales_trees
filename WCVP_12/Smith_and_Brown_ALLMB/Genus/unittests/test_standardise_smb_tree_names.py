import unittest
import re

from Smith_and_Brown_ALLMB.Genus.standardise_smb_tree_names import substitute_name_in_tree, get_binomial_from_label


class TestSubstituteNameInTreeFromRow(unittest.TestCase):
    def setUp(self):
        self.tree_string = "(Homo sapiens:0.024,(Gorilla gorilla:0.07,Pan troglodytes:0.09);"
        self.old_name = "Homo sapiens"
        self.new_name = "Homo neanderthalensis"
        self.longtree = "(Homo sapiens indicus:0.024,(Gorilla gorilla:0.07,Pan troglodytes:0.09)Homo sapines-subsp)((Homo.sapines;"

    def test_substitute_name_in_tree_from_row(self):
        expected_output = "(Homo_neanderthalensis:0.024,(Gorilla gorilla:0.07,Pan troglodytes:0.09);"
        self.assertEqual(substitute_name_in_tree(self.tree_string, self.old_name, self.new_name),
                         expected_output)

    def test_substitute_name_in_tree_from_row_no_match(self):
        no_match_string = self.tree_string
        self.assertEqual(substitute_name_in_tree(no_match_string, "No Match", self.new_name),
                         no_match_string)

    def test_substitute_name_in_tree_from_row_empty_string(self):
        self.assertEqual(substitute_name_in_tree("", self.old_name, self.new_name),
                         "")

    def test_substitute_name_in_tree_from_row_replace_with_space_in_new_name(self):
        expected_output = "(Homo_neanderthalensis_indicus:0.024,(Gorilla gorilla:0.07,Pan troglodytes:0.09);"
        self.assertEqual(substitute_name_in_tree(self.tree_string, self.old_name, "Homo neanderthalensis indicus"),
                         expected_output)

    def test_longer_examples(self):
        # Check doesn't incorrectly match
        self.assertEqual(self.longtree,substitute_name_in_tree(self.longtree, self.old_name, self.new_name))
        self.assertEqual(self.longtree, substitute_name_in_tree(self.longtree, 'Homo sapines', self.new_name))
        self.assertEqual(self.longtree, substitute_name_in_tree(self.longtree, 'Homo', self.new_name))
        self.assertEqual("(Homo sapiens indicus:0.024,(Gorilla gorilla:0.07,Pan troglodytes:0.09)Homo sapines-subsp)((Gorilla_gorilla;",
                         substitute_name_in_tree(self.longtree, 'Homo.sapines', 'Gorilla gorilla'))

    def test_problems(self):
        input = '2,Tylophora_aff._rotundifolia_Sebastian_s.n.:5.06082,V'
        new_name = 'NON_FAMILY_TIP'
        output = substitute_name_in_tree(input, 'Tylophora_aff._rotundifolia_Sebastian_s.n.', new_name)
        self.assertEqual(output, '2,NON_FAMILY_TIP:5.06082,V')

class TestGetBinomialFromLabel(unittest.TestCase):

    def test_get_binomial_from_label(self):
        # Vanilla Test Case: Both Genus and Species present, with additional indecipherable information.
        test_label_1 = "Canis_lupus_subsp."
        expected_result_1 = "Canis lupus"
        self.assertEqual(get_binomial_from_label(test_label_1), expected_result_1)

        # Test Case with different epithet indicator (var.)
        test_label_2 = "Homo_sapiens_var."
        expected_result_2 = "Homo sapiens"
        self.assertEqual(get_binomial_from_label(test_label_2), expected_result_2)

        # Test Case with no indicator and no additional information
        test_label_3 = "Felis_catus"
        expected_result_3 = "Felis catus"
        self.assertEqual(get_binomial_from_label(test_label_3), expected_result_3)

        # Test Case: Label with only genus (and no epithet), with Epithet indicator.
        test_label_4 = "Plasmodium_subsp."
        expected_result_4 = "Plasmodium"
        self.assertEqual(get_binomial_from_label(test_label_4), expected_result_4)

        # Test Case: Label have fully indecipherable information
        test_label_5 = "Anopheles_stephensi_var."
        expected_result_5 = "Anopheles stephensi"
        self.assertEqual(get_binomial_from_label(test_label_5), expected_result_5)

if __name__ == '__main__':
    unittest.main()
