import matplotlib.pyplot as plt
from ase.io.trajectory import Trajectory
from ase import Atoms
from ase.io import read
import numpy as np


def vib_analysis(model):
    '''
    Returns a list displacement of bonds/atoms in a trajectory.

    Parameters:
        model: Trajectory File

    TODO: - Resolve for periodic systems - functionally currently doesn't work


    '''
    atot = Atoms.get_chemical_symbols(self=read(model))

    traj = Trajectory(model)

    for i in range(len(atot) - 1):
        for j in range(i + 1, len(atot)):
            distances = []
            for atoms in traj:
                dist = atoms.get_distances(i, j, mic=True)
                distances.append(float(dist))
            dist_list = distances
            return dist_list


class plot_vibration_data:
    '''
    Returns a graph showing displacement of bonds/atoms in a vibration trajectory from ASE.

    Parameters:
        x_axis: length of data returned from vib_analysis()
        y_axis: data returned by vib_analysis()
        title: title for plot

    TODO: would be nice to able to plot atomic/elemental information on plot i.e which atoms are being displaced
    '''


    def __init__(self,x_axis, y_axis, title):
        self.x_axis = x_axis
        self.y_axis = y_axis
        self.title = title

    def plot_vib(self):
        plt.plot(self.x_axis, self.y_axis)
        plt.title(self.title)
        plt.show()

def vib_disps(traj, tolerance):
    '''
    Find atoms that displace during the vibration

    Parameters:
        traj: (list of atoms objects) List of atoms objects from ASE vibration .traj file
        tolerance: (float) Minimum displacement (in Angstroms) to count (Guidelines: 0.05 = main atoms, 0.01 = all atoms)

    Returns:
        Atoms and Displacements above tolerance
    '''

    max_positions, original_positions = frame_finder(traj)
    displacements = []
    for atom in range(len(original_positions)):
        move = max_positions[1].positions[atom] - max_positions[0].positions[atom]
        disp = np.linalg.norm(move)
        if disp >= tolerance:
            displacements.append([atom, disp])
    return displacements


def vib_angle(traj, vib_atom_index, ref_atom_index):
    # Only single atom vibrations supported
    traj_obj = read(traj)
    max_positions, original_positions = frame_finder(traj_obj)
    v1 = max_positions[1].positions[vib_atom_index] - max_positions[1].positions[ref_atom_index]
    v1 = v1/np.linalg.norm(v1)
    v2 = max_positions[0].positions[vib_atom_index] - max_positions[0].positions[ref_atom_index]
    v2 = v2/np.linalg.norm(v2)
    angle = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    angle_deg = angle * 180/np.pi
    return angle_deg


def characterise_vib(traj, atom_pairs):
    '''
    Characterise the vibration as (a)symmetric stretching, or bending

    Parameters:
        traj: (list of atoms objects) List of atoms objects from ASE vibration .traj file
        atom_pairs: (list of lists) list of pairs of indices of atoms forming the whole vibration

    Returns:
        Description of vibration, Full info of vibrational modes

    TODO: Add in detailed bending description (in-plane rocking/scissoring and out-of-plane wagging/twisting)
    '''
    max_positions, original_positions = frame_finder(traj)
    pair_vibs = []

    for pair in atom_pairs:
        pair_info = [pair]
        angle = vib_angle(traj, pair[0], pair[1])
        if angle > 10:
            pair_info.append('Bend')
            pair_info.append(1)
        else:
            pair_info.append('Stretch')

            if max_positions[0].get_distance(pair[0], pair[1]) - max_positions[1].get_distance(pair[0], pair[1]) > 0:
                pair_info.append(1)
            else:
                pair_info.append(0)
        pair_vibs.append(pair_info)

    if all(angle == 'Bend' for [_, angle, _] in pair_vibs):
        vib_desc = 'Bending'
    elif all(angle == 'Stretch' for [_, angle, _] in pair_vibs):
        if all(symm == 1 for [_, _, symm] in pair_vibs) or all(symm == 0 for [_, _, symm] in pair_vibs):
            vib_desc = 'Symmetric Stretching'
        else:
            vib_desc = 'Asymmetric Stretching'
    else:
        vib_desc = 'Inconclusive'

    return vib_desc, pair_vibs


def frame_finder(traj_obj):
    frames = len(traj_obj)
    original_positions = traj_obj[0]
    max_positions = [traj_obj[int(frames*1/4)], traj_obj[int(frames*3/4)]]
    return max_positions, original_positions
