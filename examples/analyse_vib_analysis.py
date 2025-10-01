'''Example of script to analyse vibrational trajectories
from a vibrational calculation using ASE'''

def test_vib_analysis():
    from carmm.analyse.vibrations import vib_analysis, plot_vibration_data, vib_disps, characterise_vib
    from ase.io import read
    from ase.visualize import view

    # Get H2O vib data calculated using ASE (EMT)
    file = 'data/H2O_vib/vib.1.traj'

    vib_data = vib_analysis(file)

    ## plot data using plot class

    x = range(len(vib_data))
    y = vib_data

    plot = plot_vibration_data(x,y,'example vib_analysis')
    plot.plot_vib()

    # vib_disp test
    traj = read('data/H2O_vib/vib.8.traj@:')
    disp_default = vib_disps(traj)
    assert disp_default == [[1, 0.352670394373909], [2, 0.35267039437390707]]

    # Displacement tolerance test
    disp_tol_test = vib_disps(traj, tolerance=0)
    assert disp_tol_test == [[0, 0.02877014449485149], [1, 0.352670394373909], [2, 0.35267039437390707]]

    # Symmetric Stretching
    traj = read('data/H2O_vib/vib.8.traj@:')
    vib_desc, pair_vibs = characterise_vib(traj, [[0, 1], [0, 2]])
    assert pair_vibs == [[[0, 1], 'Stretch', 0], [[0, 2], 'Stretch', 0]]
    assert vib_desc == 'Symmetric Stretching'

    # characterise_vib test
    # Asymetric Stretching
    traj = read('data/H2O_vib/vib.7.traj@:')
    vib_desc, pair_vibs = characterise_vib(traj, [[0, 1], [0, 2]])
    assert pair_vibs == [[[0, 1], 'Stretch', 0], [[0, 2], 'Stretch', 1]]
    assert vib_desc == 'Asymmetric Stretching'

    # Bending
    traj = read('data/H2O_vib/vib.6.traj@:')
    vib_desc, pair_vibs = characterise_vib(traj, [[0, 1], [0, 2]])
    assert pair_vibs == [[[0, 1], 'Bend', 1], [[0, 2], 'Bend', 1]]
    assert vib_desc == 'Bending'

    # Bending tolerance test
    vib_desc, pair_vibs = characterise_vib(traj, [[0, 1], [0, 2]], bending_tol=160)
    assert vib_desc == 'Symmetric Stretching'

test_vib_analysis()
