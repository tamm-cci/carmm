
def append_bulk(surface, bulk, nlayers=None, vac=10.0, output="apended.traj"):
    
    """
    Small function to append layers of bulk underneath a relaxed surface.
    Author: Hermione Beer, Oscar van Vuren

    Args:
    	surface: file or atoms containing relaxed and optimised surface geometry (.traj).
    	bulk: file or atoms containing relaxed and optimised bulk geometry (.traj). Can contain the number of layers of bulk required or be constructed (nlayers).
        nlayers: number of layers of bulk to append. If None, bulk will be appended as is
    	sep: distance (Angstrom) between two layers. If sep = 0, then the surface will sit directly on bulk with the distance between them the same as the distance between the first upper 2 layers of bulk
    	vac: vacuum to be added in the z direction
        output: name of file that will be saved (.traj) of appended bulk and surface.
    Returns:
    """
    
    from ase.io import read
    import numpy as np

    # Defining surface
    if type(surface) == str:
        surface = read(surface)
    else:
        surface=surface
    # Defining bulk
    if type(bulk) == str:
        bulk = read(bulk)
    else:
        bulk=bulk

    if nlayers is not None:
        bulk *= np.array([1,1,nlayers])
    
    # Appending them
    # Getting maximum coords of bulk in z direction
    bulk_max_z = np.max(bulk.positions[:, 2])
    # Getting minimum coords of surf in z direction
    surf_min_z = np.min(surface.positions[:, 2])
    
    # Defining vertical distance between bulk and surf
    # finding interlayer separation in the bulk
    zeds = np.sort(np.unique(bulk.positions[:,2]))
    separation = zeds[1] - zeds[0]
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
