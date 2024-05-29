import MASSiveChem.MASSiveChem as MC
import unittest

class TestFunctionalGroupDisplay(unittest.TestCase):
    def test_functional_group_display(self):
        groups_list = ['Alcohol', 'Aldehyde', 'Ketone']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(table.columns[0].field, "groups")  # Check if the first column is for functional groups
        self.assertEqual(table.columns[1].field, "images")  # Check if the second column is for images

    def test_functional_group_display_empty_list(self):
        table = MC.functional_group_display([])

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(len(table.source.data['groups']), 0)  # Check if the number of groups is 0

    def test_functional_group_display_invalid_group(self):
        groups_list = ['Alcohol', 'InvalidGroup']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)  # Check if there are two columns
        self.assertEqual(len(table.source.data['groups']), 2)  # Check if the number of groups is 1
        self.assertEqual(table.source.data['groups'][0], 'Alcohol')  # Check if the valid group is displayed

    # Additional tests
    def test_functional_group_display_multiple_valid_groups(self):
        groups_list = ['Alcohol', 'Ketone', 'Amine', 'Benzene']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 4)
        self.assertIn('Alcohol', table.source.data['groups'])
        self.assertIn('Ketone', table.source.data['groups'])
        self.assertIn('Amine', table.source.data['groups'])
        self.assertIn('Benzene', table.source.data['groups'])

    def test_functional_group_display_all_invalid_groups(self):
        groups_list = ['InvalidGroup1', 'InvalidGroup2']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 2)

    def test_functional_group_display_single_group(self):
        groups_list = ['Alcohol']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 1)
        self.assertEqual(table.source.data['groups'][0], 'Alcohol')

    def test_functional_group_display_mixed_valid_and_invalid_groups(self):
        groups_list = ['Alcohol', 'InvalidGroup', 'Benzene', 'AnotherInvalidGroup']
        table = MC.functional_group_display(groups_list)

        self.assertIsNotNone(table)
        self.assertEqual(len(table.columns), 2)
        self.assertEqual(len(table.source.data['groups']), 4)
        self.assertIn('Alcohol', table.source.data['groups'])
        self.assertIn('Benzene', table.source.data['groups'])

if __name__ == '__main__':
    unittest.main()