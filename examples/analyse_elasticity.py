def test_analyse_elasticity():
    '''
    Test the elasticity functionality
    '''
    import numpy as np
    from carmm.analyse.elasticity import read_strain_tensor_from_pkl, read_stress_from_outputs, compute_elasticity_tensor
    import os
    base = os.getcwd()
    example_path = f'{base}/data/elasticity_tensor_workflow/'
    
    # These are just checks that the file read functionality is valid
    strain_tensor = read_strain_tensor_from_pkl(example_path+'strain_tensor.npz')

    assert strain_tensor.shape == (12, 3, 3)
    assert np.isclose(np.sum(strain_tensor), 0.36) 
 
    stress_tensor = read_stress_from_outputs(path=example_path,output_file_type='.xyz')
    assert stress_tensor.shape == (12, 3, 3)

    assert np.isclose(np.sum(stress_tensor), 0.7399586028807125)
    
    # Compute elasticity tensor
    elasticity_tensor = compute_elasticity_tensor(strain_tensor, stress_tensor=stress_tensor, path=example_path, tol=1e-20)
    assert elasticity_tensor.shape == (3, 3, 3, 3)

    assert np.isclose(np.sum(elasticity_tensor), 56.758888874202626)
    assert np.isclose(elasticity_tensor[0,0,0,0], 2.87069649301488)


    # Manually save files to check functionality
    from carmm.analyse.elasticity import write_elasticity_output, write_elasticity_tensor_pickle 
    import os

    write_elasticity_output(stress_tensor, strain_tensor, elasticity_tensor, example_path)
    assert(os.path.exists(f'{example_path}/elasticity_tensor_calculation_output.txt')) # check for elasticity tensor pkl file

    write_elasticity_tensor_pickle(elasticity_tensor, example_path)

    assert (os.path.exists(f'{example_path}/elasticity_tensor.npz'))  # check for output txt file
    ########

# Run the example
test_analyse_elasticity()

def test_get_deformed_structures():

    from carmm.build.get_deformed_structures_for_elasticity_tensor import generate_deformed_strutures, create_files_and_directories
    from ase.io import read
    import os
    base = os.getcwd()

    example_path = f'{base}/data/elasticity_tensor_workflow/'
    eq_bulk = read(f'{example_path}/Co_Opt_mace_mp.traj') # reading MACE-MP optimized bulk structure.

    structure, deformations = generate_deformed_strutures(atoms_object=eq_bulk, path=example_path)
    create_files_and_directories(structure, deformations, path=example_path)

    # example of MACE-MP calculation on the deformed structures. Uncomment the lines below if you want to test it.
    # make sure MACE-MP is installed :)
    # import os
    # for defor in os.listdir(example_path):
    #     if os.path.isdir(f'{example_path}/{defor}'):
    #         atoms = read(f'{example_path}/{defor}/geometry.in')
    #         from mace.calculators import mace_mp
    #         calc = mace_mp(model='large', device='cpu', default_dtype='float64')
    #         atoms.set_calculator(calc)
    #         print(f'{defor}, {atoms.get_potential_energy()}')
    #         atoms.write(f'{example_path}/{defor}/{defor}.xyz')

    #### Assertion test for save ####

    import os
    assert(os.path.exists(f'{example_path}/defor_1/geometry.in')) # whether the deformed structures are generated
    assert (os.path.exists(f'{example_path}/strain_tensor.npz'))  # whether the strain tensor is written

test_get_deformed_structures()