def test_get_deformed_structures():

    from carmm.build.get_deformed_structures_for_elasticity_tensor import generate_deformed_strutures, create_files_and_directories
    from ase.io import read
    import os
    base = os.getcwd()

    example_path = f'{base}/data/elasticity_tensor_workflow/'
    eq_bulk = read(f'{example_path}/Co_Opt_mace_mp.traj') # reading MACE-MP optimized bulk structure.

    structure, deformations = generate_deformed_strutures(atoms_object=eq_bulk, path=example_path)
    create_files_and_directories(structure, deformations, path=example_path, copy_input_and_submission=True)

    # example of MACE-MP calculation on the deformed structures. Uncomment the lines below if you want to test it.
    # make sure MACE-MP is installed :)
    import os
    for defor in os.listdir(example_path):
        if os.path.isdir(f'{example_path}/{defor}'):
            atoms = read(f'{example_path}/{defor}/geometry.in')
            from mace.calculators import mace_mp
            calc = mace_mp(model='large', device='cpu', default_dtype='float64')
            atoms.set_calculator(calc)
            print(f'{defor}, {atoms.get_potential_energy()}')
            atoms.write(f'{example_path}/{defor}/{defor}.xyz')

    #### Assertion test for save ####

    import os
    assert(os.path.exists(f'{example_path}/defor_1/geometry.in')) # whether the deformed structures are generated
    assert (os.path.exists(f'{example_path}/strain_tensor.npz'))  # whether the strain tensor is written

