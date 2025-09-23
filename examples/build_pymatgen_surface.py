'''
This example shows how to generate a surface for rock salt MgO using the pymatgen_surface functionality
This is particularly useful when making terminations for a given facet. In this example, we will generate the 
(111) surface of rock salt MgO which exhibits two types of terminations: O-termination and Mg-termination.
'''

def test_build_pymatgen_surface():

    #### Functionality to create all surfaces ####
    from carmm.build.pymatgen_surface import generate_pymatgen_surface
    from ase.build import bulk
    bulk_model = bulk('MgO', crystalstructure='rocksalt', cubic=True, a=4.21)

    surface = generate_pymatgen_surface(bulk_model, layers=8, miller_index=(1,1,1), save=True, path='data')

    #### Assertion test for save ####
    import os
    assert(os.path.exists('data/sym_slab0_8_layer/geometry.in'))
    #########
    #
    #### Assertion tests ####
    from ase.io import read

    # ### checking slab termination
    # atoms = read('data/sym_slab0_8_layer/geometry.in')
    # z_coordinates = atoms.get_positions()[:,2]
    # import numpy as np
    # # the top and bottom atom must be oxygen as it is an oxygen terminated slab
    # assert atoms[np.argmax(z_coordinates)].symbol == 'O'
    # assert atoms[np.argmin(z_coordinates)].symbol == 'O'
    #
    # atoms = read('data/sym_slab1_8_layer/geometry.in')
    # z_coordinates = atoms.get_positions()[:, 2]
    # import numpy as np
    # # Here, the top and bottom atom must be Mg as it is a Magnesium terminated slab
    # assert atoms[np.argmax(z_coordinates)].symbol == 'Mg'
    # assert atoms[np.argmin(z_coordinates)].symbol == 'Mg'

    assert(len(atoms) == 15)

# Run the example
test_build_pymatgen_surface()
