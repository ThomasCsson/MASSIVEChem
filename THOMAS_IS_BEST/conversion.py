functional_groups_images = {
        'Alcohol': 'data/Functional groups images/Alcohol_image.png',
        'Aldehyde': 'data/Functional groups images/Aldehyde_image.png',
        'Ketone': 'data/Functional groups images/Ketone_image.png',
        'Carboxylic Acid': 'data/Functional groups images/Acid_image.png',
        'Ester': 'data/Functional groups images/Ester_image.png',
        'Ether': 'data/Functional groups images/Ether_image.png',
        'Amide': 'data/Functional groups images/Amide_image.png',
        'Amine': 'data/Functional groups images/Amine_image.png',
        'Nitrile': 'data/Functional groups images/Nitrile_image.png',
        'Chloride': 'data/Functional groups images/Halogen_image.png',
        'Bromide': 'data/Functional groups images/Bromide_image.png',
        'Fluoride': 'data/Functional groups images/Fluoride_image.png',
        'Iodide': 'data/Functional groups images/Iodide_image.png',
        'Alkene': 'data/Functional groups images/Alkene_image.png',
        'Alkyne': 'data/Functional groups images/Alkyne_image.png',
        'Imine': 'data/Functional groups images/Imine_image.png',
        'Amino acid': 'data/Functional groups images/Amino_acid_image.png',
        'Proline': 'data/Functional groups images/Proline_image.png',
        'Thiol': 'data/Functional groups images/Thiol_image.png',
        'Sulfides': 'data/Functional groups images/Sulfides_image.png',
        'Acyl Chloride': 'data/Functional groups images/Acyl_chloride_image.png',
        'Anhydride': 'data/Functional groups images/Anhydride_image.png',
        'Nitro': 'data/Functional groups images/Nitro_image.png',
        'Enamine': 'data/Functional groups images/Enamine_image.png',
        'Enamine2': 'data/Functional groups images/Enamine2_image.png',
        'Enamine3': 'data/Functional groups images/Enamine3_image.png',
        'Imide': 'data/Functional groups images/Imide_image.png',
        'Azide': 'data/Functional groups images/Azide_image.png',
        'Enol': 'data/Functional groups images/Enol_image.png',
        'Hemiacetal': 'data/Functional groups images/Hemiacetal_image.png',
        'Carbonate': 'data/Functional groups images/Carbonate_image.png',
        'Carbonate2': 'data/Functional groups images/Carbonate2_image.png',
        'Disulfide': 'data/Functional groups images/Disulfide_image.png',
        'Sulfoxide': 'data/Functional groups images/Sulfoxide_image.png',
        'Sulfone': 'data/Functional groups images/Sulfone_image.png',
        'Sulfonic acid': 'data/Functional groups images/Sulfonic_acid_image.png',
        'Thioester': 'data/Functional groups images/Thioester_image.png',
        'Phosphine': 'data/Functional groups images/Phosphine_image.png',
        'Phosphate ester': 'data/Functional groups images/Phosphate_image.png',
        'Benzene': 'data/Functional groups images/Benzene_image.png',
        'Peroxide': 'data/Functional groups images/Peroxide_image.png'
}

import base64
import pickle

def encode_image_to_base64(image_path):
    """
    Encodes an image to a base64 string.

    Args:
    image_path (str): The file path to the image.

    Returns:
    str: Base64 encoded string of the image.
    """
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')

def encode_images_in_dict(image_dict, output_file):
    """
    Encodes images in a dictionary to base64 format.

    Args:
    image_dict (dict): Dictionary with image paths as values.

    Returns:
    dict: A new dictionary with base64 encoded images as values.
    """
    base64_encoded_dict = {}
    for key, value in image_dict.items():
        try:
            with open(value, "rb") as image_file:
                base64_encoded_dict[key] = base64.b64encode(image_file.read()).decode('utf-8')
        except FileNotFoundError:
            print(f"Error: File '{value}' not found for key '{key}'.")
            base64_encoded_dict[key] = None
    
    with open(output_file, 'wb') as f:
        pickle.dump(base64_encoded_dict, f)

output_file = 'base64_encoded_images.pkl'
encode_images_in_dict(functional_groups_images, output_file)
 