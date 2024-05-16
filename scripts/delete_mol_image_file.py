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

