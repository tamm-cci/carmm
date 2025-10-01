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
    # Add an if statement to make this function also works for bcc


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
            elif abs(atoms[i].tag-atom_j.tag) > 4:
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
    :param atoms: Surface model (should be large enough, e.g.(3*3*3))
    :param a: lattice parameter
    :param lattice: crystal structure
    :param facet: tuple,e.g. (1,1,1)
    :param site: 'ontop', 'bridge'
    :return: gcn. Generalized coordination number
    """
    sum_fnn_cn = 0 # The sum of CN of fnns

    from ase.build import surface, bulk
    if lattice not in ['fcc', 'bcc']:
        print('Only support fcc and bcc for now')
        return 'gcn failed'
    

    copper = bulk('Cu', lattice, a=3.6, cubic=True)
    slab = surface(copper, facet, layers=20, vacuum=20)
    
    lastIndex = slab[-1].index
    size = len(slab)
    toplayerIndices = [atom for atom in slab if atom.tag == 1]
    toplayerSize = len(toplayerIndices)

    slab = slab.repeat((4,4,1))
    
    if site == 'ontop':
        cn_max = 12
        siteIndices = [lastIndex+5*size,]
    elif site == 'bridge':
        # Shift the indices to the centre of the surface
        siteIndices = [lastIndex+5*size, lastIndex-1+5*size]
        # Shift the indices to the centre of the surface to bulk interior
        innersiteIndices = [index-5*toplayerSize for index in siteIndices]
        cn, fnn = coord_number(slab, a=3.6, lattice=lattice, siteIndices=innersiteIndices)
        cn_max = len(fnn_set(fnn))
    
    cn, fnn = coord_number(slab, a=3.6, lattice=lattice, siteIndices=siteIndices)
    fnn_site = fnn_set(fnn)  # extracting the fnn of the site
    
    #  the CN of fnns
    cn_fnn_site, fnn_f = coord_number(slab, a=3.6, lattice=lattice, siteIndices=fnn_site)

    for n in cn_fnn_site:
        sum_fnn_cn += n   # calculating cn(j)

    gcn = sum_fnn_cn / cn_max   # dividing summation by cn_max (can do as cn_max is a constant)

    return gcn

