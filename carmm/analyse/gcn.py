def coord_number(atoms, indices=None):
    '''
    Calculate the coordination number and neighbouring atoms for any given atom or list of atoms

    Parameters:
    
    atoms: An ASE atoms object
        The object to be interrogated
    indices: 
        A list of index values for which the coordination number must be calculated
    Returns:
        - A list of coordination numbers for each interrogated species 
        - A list of neighbours for each interrogated species
    '''

    from carmm.analyse.neighbours import neighbours

    cn_list = []
    fnn_list = []

    if indices is None:
        raise ValueError("An integer value in a list form, or list of integers, is needed to evaluate the GCN")

    for i in indices:
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
    
def general_coord_number(lattice='fcc', element='Cu', a=3.6, facet=(1,1,1), site='ontop'):

    """
    GCN calculator. Only works with perfect metal slab from bulk structure for now. Accepting abitrary atoms object is not implemented for now.
    See: 
    Calle-Vallejo, F. (2023). Advanced Science, 10(20), 2207644. https://doi.org/10.1002/ADVS.202207644
    Zhao, Z., et al. (2016). Journal of Physical Chemistry C, 120(49), 28125–28130. https://doi.org/10.1021/ACS.JPCC.6B10155/
    
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

    cn, fnn = coord_number(slab, siteIndices)
    fnn_site = flatten_list_and_make_unique(fnn)  # extracting the fnn of the site
    
    #  the CN of fnns
    cn_fnn_site, fnn_f = coord_number(slab, fnn_site)

    for n in cn_fnn_site:
        sum_fnn_cn += n   # calculating cn(j)

    gcn = sum_fnn_cn / cn_max   # dividing summation by cn_max (can do as cn_max is a constant)
    return gcn

