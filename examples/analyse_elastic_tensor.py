def test_compute_elasticity_tensor():
    from carmm.analyse.calculate_elasticity_tensor import read_strain_tensor_from_pkl, read_stress_from_outputs, compute_elasticity_tensor

    example_path = 'data/elasticity_tensor_workflow'
    strain_tensor = read_strain_tensor_from_pkl(f'{example_path}/strain_tensor.pkl')
    stress_tensor = read_stress_from_outputs(path=example_path,output_file_type='.xyz')
    elasticity_tensor = compute_elasticity_tensor(strain_tensor, stress_tensor=stress_tensor, path=example_path)
    print(elasticity_tensor[0,0,0,0])


    #### Assertion test for save ####

    import os
    import numpy as np
    assert(np.isclose(elasticity_tensor[0,0,0,0], 2.9813036535, rtol=1e-8))
    assert(os.path.exists(f'{example_path}/elasticity_tensor.pkl')) # check for elasticity tensor pkl file
    assert(os.path.exists(f'{example_path}/elasticity_tensor_calculation_output.txt')) # check for output txt file
    
    ########

# Run the example
test_compute_elasticity_tensor()