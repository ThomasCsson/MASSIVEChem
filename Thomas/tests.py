import unittest
from ..THOMAS_IS_BEST.functions import data_list_generator  # Replace 'your_module' with the actual name of the module where your function is defined

class TestDataListGenerator(unittest.TestCase):
    def test_data_list_generator(self):
        # Call the function to get the generated lists
        mass, abundance, isotopes = data_list_generator()

        # Test assertions
        self.assertIsInstance(mass, list)
        self.assertIsInstance(abundance, list)
        self.assertIsInstance(isotopes, list)

        # Add more specific assertions as needed
        self.assertTrue(all(isinstance(m, float) for m in mass))
        self.assertTrue(all(0 <= ab <= 1 for ab in abundance))
        self.assertTrue(all(isinstance(iso, str) for iso in isotopes))

if __name__ == '__main__':
    unittest.main()




