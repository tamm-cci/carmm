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
    # print(stress_tensor)
    # assert stress_tensor.shape == (12, 3, 3)

    assert np.isclose(np.sum(stress_tensor), 0.7399586028807125)
    
    # Compute elasticity tensor
    elasticity_tensor = compute_elasticity_tensor(strain_tensor, stress_tensor=stress_tensor, path=example_path, tol=1e-20)
    assert elasticity_tensor.shape == (3, 3, 3, 3)
    # Troublesome assertions
    print(elasticity_tensor[0,0,0,0])
    print(np.sum(elasticity_tensor))


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


