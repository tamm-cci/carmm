def make_displaced_supercells(atoms_object, supercell_size:[1,1,1], displacement: 0.01):
    '''

    '''
    import numpy as np
    from phonopy import Phonopy
    from phonopy.interface.calculator import read_crystal_structure, write_crystal_structure
    import yaml
    from yaml import CLoader as Loader

    sup_matrix = np.diag(supercell_size)
    print(sup_matrix.shape)
    atoms_object.write('geometry_eq.in')
    unitcell, optional_structure_info = read_crystal_structure("geometry_eq.in", interface_mode='aims')
    det = np.round(np.linalg.det(sup_matrix))
    phonon = Phonopy(unitcell, supercell_matrix=sup_matrix)

    phonon.generate_displacements(distance=displacement)
    supercells = phonon.supercells_with_displacements
    phonon.save('phonopy_disp.yaml')
    stream = open("phonopy_disp.yaml", 'r')
    dictionary = yaml.load(stream, Loader)
    stream.close()
    dictionary['phonopy'].update([('calculator', 'aims'), ('configuration', {'create_displacements': '".true."', 'dim':
        f'{supercell_size[0]} {supercell_size[1]} {supercell_size[2]}', 'calculator': '"aims"'})])
    # dictionary['physical_unit'].update([('length', '"angstrom"'), ('force_constants', '"eV/angstrom^2"')])
    with open('phonopy_disp.yaml', 'w') as f:
        data = yaml.dump(dictionary, f, sort_keys=False)

    return det, supercells


def get_charges_and_moments(determinant, atoms_object):
    charges_sup = []
    moments_sup = []

    charges = atoms_object.get_initial_charges()
    moments = atoms_object.get_initial_magnetic_moments()

    for i in range(len(atoms_object)):
        for j in range(int(determinant)):
            charges_sup.append(charges[i])
            moments_sup.append(moments[i])

    return moments_sup, charges_sup


def creating_files_and_directories(supercells, charges, moments):
    '''

    :param atoms_object:
    :param charges:
    :param moments:
    :return:
    '''


    import os
    import shutil
    from phonopy.interface.calculator import write_crystal_structure
    from ase.io import read
    for ind, sup in enumerate(supercells):
        write_crystal_structure(f"geometry_{ind+1:03}.in", supercells[ind], interface_mode='aims')
        atoms = read(f"geometry_{ind+1:03}.in")
        atoms.set_initial_charges(charges)
        atoms.set_initial_magnetic_moments(moments)
        directory = f"disp-{ind+1:03}"
        parent_dir = os.getcwd()
        path_final = os.path.join(parent_dir, directory)
        if os.path.exists(path_final):
            shutil.rmtree(path_final)
        os.mkdir(path_final)
        atoms.write(path_final + '/geometry.in')
        shutil.copy(parent_dir + '/input.py', path_final + '/input.py')
        shutil.copy(parent_dir + '/submission.script', path_final + '/submission.script')
        os.chdir(parent_dir)
        os.remove(f"geometry_{ind+1:03}.in")
