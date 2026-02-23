def test_neighbours():
    '''
    Test neighbour function
    '''
    from carmm.analyse.neighbours import neighbours, first_nearest_neighbours_list

    # Build model
    from data.model_gen import get_example_slab as slab
    slab = slab(adsorbate=True)

    # Calculate neighbours
    all_neighbour_atoms, shell_list, selection = neighbours(slab, [13], 2, verbose=True)

    # Verify results
    assert(all_neighbour_atoms == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17])
    assert(shell_list == [[13], [1, 2, 4, 10, 11, 12, 14, 15, 16], [0, 3, 5, 6, 7, 8, 9, 17]])
    assert(selection[0].symbol == slab[1].symbol)
    assert(selection[0].position.all() == slab[1].position.all())

    # Calculate list for first nearest neighbours
    fnn_list = first_nearest_neighbours_list(slab, [13])
    assert fnn_list == [[1, 2, 4, 10, 11, 12, 14, 15, 16]] 

test_neighbours()
