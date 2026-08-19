def read_strain_tensor_from_pkl(pkl_file_path):
    """
    Load a strain tensor from a pickle (.pkl) file.

    Parameters
    ----------
    pkl_file_path : str
        Path to the pickle file containing the stored strain tensor.
        The file is expected to contain either:
        - a NumPy array of shape (N, 3, 3), where N is the number of deformed structures

    Returns
    -------
    strain_tensor : numpy.ndarray or object
        The strain tensor loaded from the pickle file. In most workflows,
        this is a Nx3x3 NumPy array representing the strain tensor.

    Notes
    -----
    - This function assumes the pickle file was created using the generate_deformed_structures() function with
      write_strain=True.
    """
    import numpy as np
    data = np.load(pkl_file_path)
    strain_tensor = data["strain_tensor"]
    return strain_tensor

def read_stress_from_outputs(path=None, output_file_type='.xyz'):
    """
    Read and extract the stress tensor from simulation output files.

    Parameters
    ----------
    path: str, optional
        Path pointing towards the directory that contains calculations on all the deformed structures (i.e working directory).
        Currently, the function expects the directory to have the following structure
        ---> defor_1
        ---> defor_2
        ...
        ---> defor_n
        where defor_1, defor_2,... contains the first-principles or force-field calculation outputs on the deformed
        structures.
        If None, the current working directory will be used.

    output_file_type : str, optional
        File extension specifying the type of output file to
        parse. Default is '.xyz'. Typical values include:
            - 'traj'        : ASE trajectory file containing stress information
            - 'xyz'        : XYZ file containing stress information
        The behavior depends on the implementation of the corresponding parser.

    Returns
    -------
    stress : numpy.ndarray
        A Nx3x3 stress tensor extracted from the N output files. Units depend on the
        underlying simulation code (e.g., for aims.out --> eV/Å³).

    Notes
    -----
    - For ASE trajectory files, the function will typically use
      `ase.io.read` and access stress via `atoms.get_stress()`.
    - For aims.out file, the function will manually parse the file to extract the stress tensor.
    """

    from ase.io import read
    import numpy as np
    import os

    if path is None:
        home = os.getcwd()
    else:
        home = path

    file_ext = ['.traj', '.xyz']

    if output_file_type in file_ext:
        stress_list = []
        for defor in os.listdir(home):
            if os.path.isdir(f'{home}/{defor}'):
                for file in os.listdir(f'{home}/{defor}'):
                    if file.endswith(output_file_type):
                        atoms = read(f'{home}/{defor}/{file}')
                        stress = atoms.get_stress(voigt=False)
                        stress_list.append(stress)

        stress_tensor = np.array(stress_list)

        return stress_tensor

    else:
        raise ValueError(
            f'The file extension provided is not supported. Please make sure it is one of the following [.traj, .xyz, .out]')

def compute_elasticity_tensor(strain_tensor,stress_tensor,path=None,write_elasticity_tensor=False,write_output=False,tol=1e-10):
    """
    Compute the elasticity tensor from strain and stress tensors.

    Parameters
    ----------
    strain_tensor : numpy.ndarray
        Array of shape (N, 3, 3) containing the applied strain tensors,
        where N is the number of deformation structures. Each element is a
        symmetric 3×3 strain tensor.
    stress_tensor : numpy.ndarray
        Array of shape (N, 3, 3) containing the resulting stress tensors
        corresponding to each strain in `strain_tensor`. Units depend on the underlying simulation code.
    path: str, optional
        Path pointing towards the directory that is working directory for this workflow.
        In this function, the path will be used to write the elasticity tensor as a .pkl file and an output summary as
        txt file.
        If None, the current working directory will be used.
    write_elasticity_tensor : bool, optional
        If True (default), write the computed fourth-rank elasticity tensor
        (C_{ijkl}) to disk as pickle (.pkl) file.
    write_output : bool, optional
        If True (default), print summary information such as the unique values in elasticity tensor,
         complete strain and stress tensors.
    tol : float, optional
        Value under which matrix elements are ignored. Made variable to aid CI testing.

    Returns
    -------
    elasticity_tensory : numpy.ndarray
        The fourth-rank elasticity tensor (C_{ijkl})

    Notes
    -----
    - The elasticity tensor is computed by solving the linear relation:
          σ_{ij} = C_{ijkl} ε_{kl} using the Pymatgen functionality diff_fit().
    - The function assumes small-strain elasticity.
    - The unique values are printed in output summary because the crystal symmetry causes several
      elements of the tensor to be equal. Hence, the unique values come in handy when comparing with literature values.
    """
    #import numpy as np
    from pymatgen.analysis.elasticity import diff_fit
    import os
 
    elasticity_tensor = diff_fit(strain_tensor, stress_tensor, order=2, tol=tol)

    if path is None:
        path = os.getcwd()

    if write_output:
        write_elasticity_output(stress_tensor,strain_tensor,elasticity_tensor,path)

    if write_elasticity_tensor:
        write_elasticity_tensor_pickle(elasticity_tensor, path)

    return elasticity_tensor[0]

def write_elasticity_output(stress_tensor, strain_tensor, elasticity_tensor, path):
    """
    Writes an output file (.txt format) which includes summary information such as the complete strain and stress
    tensors and the unique values of elasticity tensor in the units of Pa and eV/Å³.

    Parameters
    ----------
    strain_tensor : numpy.ndarray
        Array of shape (N, 3, 3) containing the applied strain tensors,
        where N is the number of deformation structures. Each element is a
        symmetric 3×3 strain tensor.
    stress_tensor : numpy.ndarray
        Array of shape (N, 3, 3) containing the resulting stress tensors
        corresponding to each strain in `strain_tensor`. Units depend on the underlying simulation code.
    elasticity_tensor : numpy.ndarray
        The fourth-rank elasticity tensor (C_{ijkl})
    path : string
        Directory for storing files

    Returns
    -------
    None
        This function writes output text to disk.

    """

    import numpy as np
    try:
        f = open(f'{path}/elasticity_tensor_calculation_output.txt', 'x')
        f.close()
    except:
        print('The file already exists. Overwriting the elasticity_tensor_calculation_output.txt file...')

    f = open(f'{path}/elasticity_tensor_calculation_output.txt', 'w')
    f.write('Unique values in elasticity tensor in units of Pascals (Pa).\n')
    f.write(str(np.round(np.unique(np.array(elasticity_tensor[0]), return_index=True)[0], 3) * 160.2176621) + '\n')
    f.write('Unique values in elasticity tensor in units of eV/Å³.\n')
    f.write(str(np.round(np.unique(np.array(elasticity_tensor[0]), return_index=True)[0], 3)) + '\n')
    f.write('Elasticity tensor in full. The above values show the unique values in the tensor.\n')
    f.write(str(np.array(elasticity_tensor[0])) + '\n')
    f.write('Stress tensor for all deformations in full. Shape is Nx3x3 where N is the number of deformations '
            'used to compute the elasticity tensor.\n')
    f.write(str(stress_tensor) + '\n')
    f.write(str(stress_tensor.shape) + '\n')
    f.write('Strain tensor for all deformations in full. Shape is Nx3x3 where N is the number of deformations '
            'used to compute the elasticity tensor.\n')
    f.write(str(strain_tensor) + '\n')
    f.write(str(strain_tensor.shape) + '\n')
    f.close()

def write_elasticity_tensor_pickle(elasticity_tensor, path):

    """
    Writes an output file (.pkl format) which contains the elasticity tensor

    Parameters
    ----------
    elasticity_tensor : numpy.ndarray
        The fourth-rank elasticity tensor (C_{ijkl})
    path: string
        The directory where the file should be saved

    Returns
    -------
    None
        This function writes output to disk.
    """

    import numpy as np

    np.savez(
        f"{path}/elasticity_tensor.npz",
        elasticity_tensor=np.array(elasticity_tensor[0]))






