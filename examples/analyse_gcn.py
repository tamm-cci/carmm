def test_gcn():
    # This tests the GCN calculator.
    from carmm.analyse.gcn import generalised_coordination_number

    # Build model. This has a CO2 adsorbate
    from data.model_gen import get_example_slab
    slab = get_example_slab(adsorbate=False)

    # The GCN here would be 7.5 for a real system, but the slab is only 2 layers thick
    gcn = generalised_coordination_number(slab, [1], 12)
    assert gcn == 6.75    

    # To prove this, we make a thicker slab
    slab_thicker = get_example_slab(adsorbate=False, thickness=3)
    gcn = generalised_coordination_number(slab_thicker, [1], 12)
    assert gcn == 7.5

    # A note these calculations below are not physically valid, except for the "FCC" case.
    # The rest are just there to completely test the code functionality.   
    for lattice in ['sc', 'bcc', 'fcc', 'hcp']:
        gcn = generalised_coordination_number(slab_thicker, [1], lattice)
        if lattice == 'sc':
            assert gcn == 15.0
        elif lattice == 'bcc':
            assert gcn == 11.25
        else:
            # FCC and HCP
            assert gcn == 7.5

test_gcn()
