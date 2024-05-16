import os

def delete_mol_image_file():

    #---------------------------------------------------------------------------------------------#
    '''
    delete_mol_image_file()
    
    Input: None
    
    Output: Deletes the creates file to store the image : molecule_image.png
    '''
    #---------------------------------------------------------------------------------------------#

    # file path to delete
    filepath = 'molecule_image.png'

    #checks if the file exists
    if os.path.exists(filepath):

        #if yes, deletes the file
        os.remove(filepath)

        print(f"File '{filepath}' has been deleted.")
    
    else:

        print(f"File '{filepath}' does not exist.")

import unittest
from unittest.mock import patch

class TestDeleteMolImageFile(unittest.TestCase):
    @patch('builtins.print')
    @patch('os.path.exists')
    @patch('os.remove')
    def test_delete_mol_image_file(self, mock_remove, mock_exists, mock_print):
        # Test when the file exists
        mock_exists.return_value = True
        delete_mol_image_file()
        mock_remove.assert_called_once_with('molecule_image.png')
        mock_print.assert_called_once_with("File 'molecule_image.png' has been deleted.")

    @patch('builtins.print')
    @patch('os.path.exists')
    @patch('os.remove')
    def test_delete_mol_image_file_nonexistent(self, mock_remove, mock_exists, mock_print):
        # Test when the file does not exist
        mock_exists.return_value = False
        delete_mol_image_file()
        mock_remove.assert_not_called()
        mock_print.assert_called_once_with("File 'molecule_image.png' does not exist.")

if __name__ == '__main__':
    unittest.main()
