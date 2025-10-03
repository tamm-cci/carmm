
def append_bulk(input_surface, input_bulk, nlayers=None, vac=10.0, output="apended.traj"):
    
    """
    Small function to append layers of bulk underneath a relaxed surface.
    Allows for the construction of very large slab models for cutting QM/MM clusters
    without having to relax such a large model.
    Author: Hermione Beer, Oscar van Vuren

    Args:
    	surface: file or atoms containing relaxed and optimised surface geometry (.traj).
    	bulk: file or atoms containing relaxed and optimised bulk geometry (.traj). Can contain the number of layers of bulk required or be constructed (nlayers).
        NOTE: THE REPEAT UNIT OF BULK IN A SLAB MAY NOT BE THE BULK UNIT CELL. ENSURE USE OF CORRECT GEOMETRY FOR BULK
        nlayers: number of layers of bulk to append. If None, bulk will be appended as is
    	sep: distance (Angstrom) between two layers. If sep = 0, then the surface will sit directly on bulk with the distance between them the same as the distance between the first upper 2 layers of bulk
    	vac: vacuum to be added in the z direction
        output: name of file that will be saved (.traj) of appended bulk and surface.
    Returns:
    """
    
    from ase.io import read
    import numpy as np

    # Defining surface
    if type(input_surface) == str:
        surface = read(input_surface)
    else:
        surface=input_surface
    # Defining bulk
    if type(input_bulk) == str:
        bulk = read(input_bulk)
    else:
        bulk=input_bulk

    if nlayers is not None:
        layers_mult = np.array([1,1,nlayers])
        print(layers_mult)
        bulk = bulk * layers_mult
    
    # Appending them
    # Getting maximum coords of bulk in z direction
    bulk_max_z = np.max(bulk.positions[:, 2])
    # Getting minimum coords of surf in z direction
    surf_min_z = np.min(surface.positions[:, 2])
    
    # Defining vertical distance between bulk and surf
    # finding interlayer separation in the bulk
    zeds = np.sort(np.unique(bulk.positions[:,2]))
    xyzero = [bulk.positions[0,0], bulk.positions[0,1]]
    # Ensuring that we are only counting the separation in the z dimension
    for ind in range(len(zeds)):
        if bulk.positions[ind,0] == xyzero[0] and bulk.positions[ind,1] == xyzero[1]:
            separation = separation = zeds[ind] - zeds[0]
            break
        else:
            continue

    # Constructing shift using interlayer distance
    z_shift = bulk_max_z - surf_min_z + separation

    surface.positions[:, 2] += z_shift # Apply transition to surf
    appended_material = bulk + surface # Combining bulk and surf into new object
    
    # Resizing cell and centering it
    final_cell = appended_material.get_cell()
    final_cell[2, 2] = np.max(appended_material.positions[:, 2]) - np.min(appended_material.positions[:, 2]) + vac
    appended_material.set_cell(final_cell)
    appended_material.center(axis=2)
    
    # Saving it
    if output is not None:
        appended_material.write(output)

    return appended_material
