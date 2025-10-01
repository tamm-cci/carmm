def coord_number(atoms, a=3.6, lattice='fcc', siteIndices=None):
    # The list that stores coordination number for each atom
    cn_list = []
    # The list that stores the indices of first nearst neighbours for each atom
    # This is a list of fnn index lists
    fnn_list = []

    if lattice == 'fcc':
        bond = round(a / 2 ** 0.5, 3)
    elif lattice == 'bcc':
        bond = round(a * (3 ** 0.5) / 2, 3)


    for i in siteIndices:
        cn = 0
        # List for first nearest neighbours of atom i
        fnn = []
        # Counting coordination number for atom i
        for atom_j in atoms:
            j = atom_j.index
            # Skip the iteration if we are considering the distance
            # between atom i and itself
            if i == j:
                continue
            elif abs(atoms[i].z-atom_j.z) > 4:
                continue
            # Check if atom i and atom j are first nearst neighbours
            dist = atoms.get_distance(i, j, mic=True)
            if dist <= bond+0.001 and dist >= bond-0.001:
                cn += 1
                fnn.append(j)

        # Append coordination number and first nearest neighbours to the lists every time the second j loop finishes.
        cn_list.append(cn)
        fnn_list.append(fnn)

    return cn_list, fnn_list

def fnn_set(fnnLists):
    """
    :param fnnLists
    :return: Set of all fnns
    """
    allfnn = []
    for fnn in fnnLists:
        allfnn += fnn
    return set(allfnn)
    
def general_coord_number(lattice='fcc', facet=(1,1,1), site='ontop'):

    """
    GCN calculator.
    See: 
    Calle-Vallejo, F. (2023). Advanced Science, 10(20), 2207644. https://doi.org/10.1002/ADVS.202207644
    Zhao, Z., et al. (2016). Journal of Physical Chemistry C, 120(49), 28125–28130. https://doi.org/10.1021/ACS.JPCC.6B10155/
    :param atoms: Surface model (should be large enough, e.g.(3*3*3))
    :param a: lattice parameter
    :param lattice: crystal structure
    :param facet: tuple,e.g. (1,1,1)
    :param site: 'ontop', 'bridge'
    :return: gcn. Generalized coordination number
    """
    from ase.build import surface, bulk
    from ase.visualize import view
    import numpy as np

    sum_fnn_cn = 0 # The sum of CN of fnns
    a =3.6
    
    if lattice == 'fcc':
        bond = round(a / 2 ** 0.5, 3)
    elif lattice == 'bcc':
        bond = round(a * (3 ** 0.5) / 2, 3)
    
    if lattice not in ['fcc', 'bcc']:
        print('Only support fcc and bcc for now')
        return 'gcn failed'
    

    copper = bulk('Cu', lattice, a=a, cubic=True)
    slab = surface(copper, facet, layers=20, vacuum=20)
     
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
    cn, fnn = coord_number(slab, a=a, lattice=lattice, siteIndices=innersiteIndices)
    cn_max = len(fnn_set(fnn))
    print('CN-max ', cn_max)

    cn, fnn = coord_number(slab, a=a, lattice=lattice, siteIndices=siteIndices)
    fnn_site = fnn_set(fnn)  # extracting the fnn of the site
    
    #  the CN of fnns
    cn_fnn_site, fnn_f = coord_number(slab, a=a, lattice=lattice, siteIndices=fnn_site)

    for n in cn_fnn_site:
        sum_fnn_cn += n   # calculating cn(j)

    gcn = sum_fnn_cn / cn_max   # dividing summation by cn_max (can do as cn_max is a constant)
    # view(slab)
    return gcn

