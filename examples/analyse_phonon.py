# import sys
# if sys.version_info.major==3 and sys.version_info.minor!=8:
def test_phonon_analysis():
    '''
    After creating the displaced geometries in seperate folders, you will have to run first-principles calculation
    in FHI-aims in each of the separate directory to generate a aims.out file in each disp-00n folder
    The post-process of phonon calculations using phonon_analysis.py script in analyse folder depends on the presence
    of aims.out file in each of these folders

    '''

    from ase.build import bulk
    from carmm.analyse.phonon_analysis import get_band_conf, get_thermal_conf, generate_phonon_data, phonon_data_to_csv
    from ase.io import read
    import os

    example_path = 'data/phonon_workflow'

    atoms = read(f'{example_path}/geometry_eq.in')
    get_band_conf(atoms, path=example_path)
    get_thermal_conf(path=example_path)

    assert (os.path.exists(f'{example_path}/band.conf'))
    assert (os.path.exists(f'{example_path}/thermal.conf'))

    generate_phonon_data(path=example_path)
    phonon_data_to_csv(band_data=True, path=example_path)

    generated_files = ['FORCE_SETS','band.yaml','thermal_properties.yaml','band_data.csv','thermal_data.csv']

    for file in generated_files:
        assert (os.path.exists(f'{example_path}/{file}'))

from build_get_displaced_structures_phonon import test_generate_displaced_structures
test_generate_displaced_structures()
test_phonon_analysis()
