def coord_number(atoms, site_indices=None):
    '''
    Calculate the coordination number and neighbouring atoms for any given atom or list of atoms

    Parameters:
    
    atoms: An ASE atoms object
        The object to be interrogated
    site_indices: 
        A list of index values for which the coordination number must be calculated
    Returns:
        - A list of coordination numbers for each interrogated species 
        - A list of neighbours for each interrogated species
    '''

    from carmm.analyse.neighbours import neighbours

    cn_list = []
    fnn_list = []

    if site_indices is None:
        raise ValueError("An integer value in a list form, or list of integers, is needed to evaluate the GCN")

    for i in site_indices:
        all_neighbour_atoms, shell_list, selection = neighbours(atoms, [i], 1, verbose=False)
        cn_list.append(len(shell_list[1]))
        fnn_list.append(shell_list[1])

    return cn_list, fnn_list

def flatten_list_and_make_unique(list_of_lists):
    """
    This converts a list of lists into a flatten 1D list of integers with unique values

    Parameters:

    list_of_lists: List of integer lists 
        Contains the first-nearest neighbours for each atom
    Returns:
        - A flattened and unique 1D list of all possible second neighbours (neighbours of first neighbours)    
 
    """
    
    all_values = []
    for list in list_of_lists:
        all_values += list
    
    # Demonstration how to do the above in one line - saving for referece
    # all_values = [value for list in list_of_lists for value in list]

    # Make the list unique
    return set(all_values)
    
def setup_metal_slab_for_general_coordination_number_calculation(lattice='fcc', element='Cu', a=3.6, facet=(1,1,1), site='ontop'):

    """
       


    :param lattice: crystal structure
    :param facet: tuple, e.g. (1,1,1)
    :param site: 'ontop', 'bridge'
    :return: gcn. Generalized coordination number
    """
    from ase.build import surface, bulk
    import numpy as np

    sum_fnn_cn = 0 # The sum of CN of fnn
    
    if lattice == 'fcc':
        bond = round(a / 2 ** 0.5, 3)
    elif lattice == 'bcc':
        bond = round(a * (3 ** 0.5) / 2, 3)
    
    if lattice not in ['fcc', 'bcc']:
        print('Only support fcc and bcc for now')
        raise KeyError('Generalised coordination number is only supported for FCC and BCC lattice')

    metalbulk = bulk(element, lattice, a=a, cubic=True)
    slab = surface(metalbulk, facet, layers=20, vacuum=20)
     
    size = len(slab)
    maxZ = np.max([atom.z for atom in slab])

    lastIndex = 0 # Surface atom at the top right corner
    toplayerIndices = [atom.index for atom in slab if atom.z == maxZ]
    toplayerSize = len(toplayerIndices)
    
    x=0
    y=0
    for index in toplayerIndices:
        if slab[index].x > x and slab[index].y > y:
            x = slab[index].x
            y = slab[index].y
            lastIndex = index
    
    # FCC cn_max = 12
    
    if site == 'ontop':
        siteIndices = [lastIndex+5*size,]
    elif site == 'bridge':
        # Shift the indices to the centre of the surface
        topfnn = [atom.index for atom in slab 
                  if atom.z == maxZ and 
                  slab.get_distance(atom.index, lastIndex) <= bond+0.001 
                  and slab.get_distance(atom.index, lastIndex) >= bond-0.001]
        
        siteIndices = [lastIndex+5*size, topfnn[0]+5*size]

    slab = slab.repeat((4,4,1))

    print('last index ', lastIndex, 'site indices', siteIndices)

    # Shift the indices to the centre of the surface to bulk interior
    innersiteIndices = [siteIndices[0]-9*toplayerSize]
    cn, fnn = coord_number(slab, innersiteIndices)
    cn_max = len(flatten_list_and_make_unique(fnn))
    print('CN-max ', cn_max)

    return generalised_coordination_number(slab, siteIndices, cn_max)


def generalised_coordination_number(slab, site_index, cn_max): 

    """
    Calculator for the generalised coordination number of an atom on a surface. 
    Can work for _any_ system but one must know the bulk coordination number for a given species.    
 
    See:
      Calle-Vallejo, F. (2023). Advanced Science, 10(20), 2207644. https://doi.org/10.1002/ADVS.202207644
      Zhao, Z., et al. (2016). Journal of Physical Chemistry C, 120(49), 28125–28130. https://doi.org/10.1021/ACS.JPCC.6B10155/

    Parameters:
    
    slab: ASE atoms object
        Representation of the surface for consideration
    site_index: Integer
        The site of interest for which GCN should be calculated
    cn_max: Integer
        The maximum coordination number for the species being considered
    TODO: Could this also be considered a string, with values stored for common entries (FCC, BCC, etc.)?

    Returns:
      - Float value of the generalised coordination number

    """

    # Calculate first nearest neighbours, and their coordination
    _, first_nearest_neighbours = coord_number(slab, site_index)
    fnn_flattened = flatten_list_and_make_unique(first_nearest_neighbours)
    coordination_of_fnn, _ = coord_number(slab, fnn_flattened)

    # Dividing summation by cn_max (can do as cn_max is a constant)	
    gcn = sum(coordination_of_fnn) / cn_max  
    
    return gcn

