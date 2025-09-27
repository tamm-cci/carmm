def get_band_conf(atoms):
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
        # elif sp == ',':
        #     band[-1] = band[-1] + ','

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
    # f.write(f'TPROP =.TRUE.\n'
    #         f'MESH = 16 16 16\n')

def get_thermal_conf(atoms):
    f = open('thermal.conf', 'w')
    f.write(f'TPROP =.TRUE.\n'
            f'MESH = 16 16 16\n')
    f.close()

def generate_phonon_data():
    import os
    os.system("phonopy -f disp-?/aims.out")
    os.system("phonopy -p -s band.conf")
    os.system("phonopy -p -s thermal.conf")

def phonon_data_to_csv(band_data=False, thermal_data=True, band_file='band.yaml', thermal_file='thermal_properties.yaml'):
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
            raise Exception('The file \'thermal_properties.yaml\' not found. Make sure you are providiing the correct file name to the '
                  'thermal_file parameter in the function')
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
            raise Exception('The file \'band.yaml\' not found. Make sure you are providiing the correct file name to the '
                  'band_file parameter in the function')
