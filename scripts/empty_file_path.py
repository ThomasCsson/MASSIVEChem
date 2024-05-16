import os

def empty_file_path(search_directory='.'):

     #---------------------------------------------------------------------------------------------#
    '''
    empty_file_path()
    
    Input: search_directory, which specifies to check in the current directory
    
    Output: - creates a new file called molecule_image.png to store an image later
            - returns the relative path to the file
    '''
    #---------------------------------------------------------------------------------------------#

    # name of the file name to create
    filename = 'molecule_image.png'

    #finds the current directory
    current_directory = os.getcwd()

    #creates a file path to the current directory
    filepath = os.path.join(current_directory, filename)

    #checks if the file already exists
    if not os.path.exists(filepath):

        #if no, creates the path
        with open(filepath, 'a'):
            pass
    else:

        #else pass
        pass
    
    
    for root, dirs, files in os.walk(search_directory):

        #checks all the file names in the directory
        if filename in files:

            #return file path
            return os.path.join(root, filename)
    
    return None


print(empty_file_path())