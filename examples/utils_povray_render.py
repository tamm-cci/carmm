"""
Test script for povray_render and atom_sub functions
Usage example shown first, other tests shown after.
"""


def test_povray_render():
    from carmm.utils.povray_render import povray_render, atom_sub
    from ase.io import read
    import numpy as np
    from ase.data import covalent_radii, atomic_numbers
    from ase.data.colors import jmol_colors, cpk_colors

    # Main example script
    atoms = read('data/NH3-H3O_traj/nh3-h3o.traj')

    # An alternative to substituting atoms: changing the colour of specific elements (use any of the below methods)
    # Use "'colors': None," for all default colours
    color_dict = {
        'H': np.array([255, 255, 255])/255,  # 1. RGB colours converted to between 0 and 1
        'N': jmol_colors[atomic_numbers['N']],  # 2. Use default jmol_colors or cpk_colors in ASE
        # 3. Don't specify atom color: function will automatically use default jmol_colors in ASE
    }

    # Atoms can be made greyscale by converting their RGB colors using the standard RGB luminosity formula:
    # newR = newG = newB = 0.2126R + 0.7152G + 0.0722B
    # e.g. Oxygen jmol RGB = [1, 0.051, 0.051] -> Greyscale oxygen RGB = [0.2528, 0.2528, 0.2528]

    radius_scale = 0.8
    radius_list = []
    for atomic_number in atoms.get_atomic_numbers():
        if atomic_number == 8:
            radius_list.append(radius_scale * 1)
        else:
            radius_list.append(radius_scale * covalent_radii[atomic_number])

    generic_projection_settings = {
        'rotation': '90x,80y,90z',
        'radii': radius_list,
        'colors': color_dict,
    }

    povray_settings = {
        'camera_type': 'orthographic angle 7',
        'camera_dist': 100,
        'celllinewidth': 0.07,
        'textures': ['jmol'] * len(atoms),
        'transmittances': np.linspace(0, 1, len(atoms)),  # Opaque atom = 0, Transparent atom = 1
    }

    povray_render(atoms, generic_projection_settings=generic_projection_settings, povray_settings=povray_settings)

    # Tests

    # Custom settings

    gen_proj_sett, povray_sett = povray_render(atoms,
                                               generic_projection_settings=generic_projection_settings,
                                               povray_settings=povray_settings)


    assert gen_proj_sett['rotation'] == '90x,80y,90z'
    assert gen_proj_sett['radii'] == radius_list
    np.testing.assert_almost_equal(gen_proj_sett['colors']['H'], np.array([1, 1, 1]), decimal=3)
    np.testing.assert_almost_equal(gen_proj_sett['colors']['N'], np.array([0.188, 0.314, 0.973]), decimal=3)
    np.testing.assert_almost_equal(gen_proj_sett['colors']['O'], np.array([1, 0.051, 0.051]), decimal=3)

    assert povray_sett['camera_type'] == 'orthographic angle 7'
    assert povray_sett['camera_dist'] == 100
    assert povray_sett['celllinewidth'] == 0.07
    assert np.all(povray_sett['textures'] == ['jmol'] * len(atoms))
    assert povray_sett['display'] is False
    assert np.all(np.isclose(povray_sett['transmittances'],
                             [0, 0.14285714, 0.28571429, 0.42857143, 0.57142857, 0.71428571, 0.85714286, 1],
                             atol=0.001))

    # Default settings

    gen_proj_sett, povray_sett = povray_render(atoms)

    assert np.all(list(gen_proj_sett.keys()) == ['rotation', 'radii', 'colors'])
    assert gen_proj_sett['rotation'] == '0x,0y,0z'
    assert gen_proj_sett['radii'] == 1.0
    true_colors = [[1, 0.051, 0.051], [1, 1, 1], [1, 1, 1], [1, 1, 1],
                   [0.188, 0.314, 0.973], [1, 1, 1], [1, 1, 1], [1, 1, 1]]
    for i in np.arange(len(gen_proj_sett['colors'])):
        assert np.all(gen_proj_sett['colors'][i] == true_colors[i])

    assert np.all(list(povray_sett.keys()) == ['camera_type', 'camera_dist', 'transmittances', 'display'])
    assert povray_sett['camera_type'] == 'orthographic angle 5'
    assert povray_sett['camera_dist'] == 50
    assert np.all(povray_sett['transmittances'] == [0] * len(atoms))
    assert povray_sett['display'] == False

    # Warnings

    povray_render(atoms, povray_settings={'camera_type': 'orthographic'})
    povray_render(atoms, povray_settings={'camera_type': 'perspective'})
    povray_render(atoms, povray_settings={'camera_type': 'ultra_wide_angle'})

    # Atom substitution function, useful for easily changing element colours
    subbed_atoms = atom_sub(atoms, atom_subs=[['N', 'C'], ['O', 'N']])
    assert subbed_atoms.symbols[4] == 'C'
    assert subbed_atoms.symbols[0] == 'N'
    assert len(subbed_atoms) == 8


test_povray_render()
