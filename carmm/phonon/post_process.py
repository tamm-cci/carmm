def get_band_conf(atoms):
    """
    Generate the band structure configuration (band.conf file) for a given atomic structure.
    The band.conf file is necessary to generate a phonon bandstructure once the force set are collected.

    Parameters
    ----------
    atoms : ase.Atoms
        structure of the unit cell as ase.Atoms object to generate the high symmetry points in the reciprocal space.

    Returns
    -------
    None
    Generates a band.conf file in the working directory that will be required for further post-processing.

    Notes
    -----
    - This function serves as a helper to prepare configuration file to generate the phonon bandstructure.
    - The exact k-path convention may depend on the crystal system

    """


    lat1 = atoms.cell.bandpath()
    band_label = []
    band = []
    for sp in lat1.path:
        if sp == 'G':
            band_label.append('$\\Gamma$')
        elif sp.isnumeric():
            band_label[-1] = band_label[-1] + sp
        else:
            band_label.append(sp)
    for sp in band_label:
        if sp != ',':
            if sp == '$\\Gamma$':
                band.append('0 0 0')
            else:
                cor1 = lat1.special_points[sp][0]
                cor2 = lat1.special_points[sp][1]
                cor3 = lat1.special_points[sp][2]
                band.append(f'{cor1} {cor2} {cor3}')

    for i in band_label:
        if i == ',':
            band_label.pop(band_label.index(i))

    band_label_final = ' '.join(band_label)
    print(band_label_final)
    band_final = '    '.join(band)
    print(band_final)

    try:
        f = open('band.conf', 'x')
        f.close()
    except:
        print()

    f = open('band.conf', 'w')
    f.write(f'BAND = {band_final}\n'
            f'BAND_LABELS = {band_label_final}\n'
            f'BAND_POINTS = 101\n')
    f.close()


def get_thermal_conf():
    """
    Generate the thermal properties configuration (thernal.conf file)
    The thermal.conf file is necessary to generate a thermal properties data once the force set are collected and the
    phonon bandstructure is generated.

    Parameters
    ----------
    None

    Returns
    -------
    None
    Generates a thermal.conf file in the working directory that will be required for further post-processing.

    Notes
    -----
    - This function serves as a helper to prepare configuration file to generate the thermal property data like
    energy, entropy, Helmholtz free energy and heat capacity.
    - NOTE: currently the temperature range is from 0 to 1000 K. For a majority of systems, the range is sufficient to
    understand the thermal properities. Future development may allow for custom range.

    """
    f = open('thermal.conf', 'w')
    f.write(f'TPROP =.TRUE.\n'
            f'MESH = 16 16 16\n')
    f.close()


def generate_phonon_data(bandstructure=True, thermal_properties=True):
    """
    Generate phonon-related data such as phonon band structures and
    thermodynamic properties.

    Parameters
    ----------
    bandstructure : bool, optional
        If True (default), compute the phonon band structure along a chosen
        high-symmetry path in the Brillouin zone using the information available in the band.conf file obtained from
        get_band_conf function.

    thermal_properties : bool, optional
        If True (default), compute thermodynamic quantities derived from the phonon density of states (DOS), such as
        free energy, entropy, and heat capacity as a function of temperature. This needs the thermal.conf file
        generated using the get_thermal_conf function.

    Returns
    -------
    None.

    Runs the terminal commands and generates the following files:
    FORCE_SETS -- force information obtained using the phonopy API. This is needed for all post-processing of data
    band.pdf -- phonon band structure
    band.yaml -- phonon band structure data as yaml file. use the phonon_data_to_csv to obtain the data as a csv file
    thermal.pdf -- thermal properties plotted as a function of temperature
    thermal_properites.yaml - thermal properties data as yaml file. use the phonon_data_to_csv to obtain the data as a
    csv file


    Notes
    -----
    To obtain the thermal properities, it is necessary to keep the paramter bandstructure=True, because thermal
    properties are derived from phonon density of states (DOS). In that case, if thermal_properties=True and bandstructure=False,
    the code will run the bandstructure command automatically to avoid any errors.
    """
    import os
    os.system("phonopy -f disp-*/aims.out")
    if bandstructure:
        os.system("phonopy -p -s band.conf")
        if thermal_properties:
            os.system("phonopy -p -s thermal.conf")
    else:
        if not bandstructure:
            if thermal_properties:
                os.system("phonopy -p -s band.conf")
                os.system("phonopy -p -s thermal.conf")


def phonon_data_to_csv(band_data=False, thermal_data=True,
                       band_file='band.yaml', thermal_file='thermal_properties.yaml'):
    """
    Convert phonon calculation results from YAML files into CSV format.

    Parameters
    ----------
    band_data : bool, optional
        If True, parse phonon band structure data from `band_file` (needs to be .yaml) and
        export it to a CSV file. Default is False.
    thermal_data : bool, optional
        If True, parse phonon thermal property data from `thermal_file` (needs to be .yaml) and
        export it to a CSV file. Default is True.
    band_file : str, optional
        Path to the YAML file containing phonon band structure data.
        Default is 'band.yaml'.
    thermal_file : str, optional
        Path to the YAML file containing thermal property data.
        Default is 'thermal_properties.yaml'.

    Returns
    -------
    None
        The function writes CSV files to disk and does not return any value.

    Notes
    -----
    - The band structure YAML file (e.g., produced by Phonopy) typically
      contains q-points and frequencies
    - The thermal properties YAML file contains temperature-dependent
      quantities such as free energy, entropy, and heat capacity.
    - The resulting CSV files are formatted for easy plotting or further
      numerical analysis with standard tools.
    """

    import yaml
    from yaml import load
    from yaml import CLoader as Loader
    import os
    import csv
    if thermal_data:
        if os.path.exists(thermal_file):
            property_list = ['Energy (kJ/mol)', 'Entropy (J/K/mol)', 'Helmholtz free energy (kJ/mol)', 'Heat Capacity (J/K/mol)']
            header = ['Temperature (K)']+property_list
            csv_file = f'thermal_data.csv'
            with open(csv_file, mode='w', newline='') as file:
                writer = csv.writer(file)
                writer.writerow(header)
                stream = open(thermal_file, 'r')
                dictionary = yaml.load(stream, Loader)
                for data in dictionary['thermal_properties']:
                    writer.writerow([data['temperature'], data['energy'], data['entropy'], data['free_energy'],
                                     data['heat_capacity']])
                stream.close()
                print('Thermal data written to file thermal_data.csv')
        else:
            raise Exception('The file \'thermal_properties.yaml\' not found. Make sure you are providing the correct '
                            'file name to the thermal_file parameter in the function')
    if band_data:
        if os.path.exists(band_file):
            csv_file = f'band_data.csv'
            with open(csv_file, mode='w', newline='') as file:
                stream = open(band_file, 'r')
                dictionary = yaml.load(stream, Loader)
                freqencies_header = [f'frequency {i+1}' for i in range(len(dictionary['phonon'][0]['band']))]
                header = ['q-position[0]', 'q-position[1]', 'q-position[2]', 'Distance']+freqencies_header
                writer = csv.writer(file)
                writer.writerow(header)
                for key in dictionary['phonon']:
                    row = [key['q-position'][i] for i in range(3)]+[key['distance']]+[j['frequency']for j in key['band']]
                    writer.writerow(row)
                stream.close()
                print('Phonon band structure data written to file band_data.csv')
        else:
            raise Exception('The file \'band.yaml\' not found. Make sure you are providing the correct file name to '
                  'the band_file parameter in the function')
