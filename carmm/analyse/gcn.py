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
    cn_max: String or Integer
        The crystal structure of the lattice (so value taken from stored information) or 
        the maximum coordination number for the species being considered.

        Accepted string values are: SC, BCC, FCC and HCP.

    Returns:
      - Float value of the generalised coordination number

    """

    from carmm.analyse.neighbours import first_nearest_neighbours_list

    # Manage the definition of cn_max; if string convert to int
    if isinstance(cn_max, str):
        if cn_max.lower() == "sc":
            cn_max = 6
        elif cn_max.lower() == "bcc":
            cn_max = 8
        elif cn_max.lower() == "fcc" or cn_max.lower() == "hcp":
            cn_max = 12
        else:
            print(cn_max, "is not recognised. Please just give the bulk coordination number as an integer.")
            raise ValueError('Tabulated coordination numers only exist for simple cubic, BCC, FCC, and HCP.') 

    # Calculate first nearest neighbours, and their coordination
    first_nearest_neighbours = first_nearest_neighbours_list(slab, site_index)
    fnn_flattened = flatten_list_and_make_unique(first_nearest_neighbours)
    list_of_first_neighbours_of_first_neighbours = first_nearest_neighbours_list(slab, fnn_flattened)
    cn_first_nearest_neighbours = [len(neighbours) for neighbours in list_of_first_neighbours_of_first_neighbours]
    sum_cn_fnn = sum(cn_first_nearest_neighbours)

    # Dividing summation by cn_max (can do as cn_max is a constant)	
    return sum_cn_fnn / cn_max  
 
def flatten_list_and_make_unique(list_of_lists):
    """
    This converts a list of lists into a flatten 1D list of integers with unique values
    TODO: Generalise and move to utils

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
    
