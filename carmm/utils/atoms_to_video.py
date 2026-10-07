from ase.io import read
from ase.visualize import view
from ase.io.trajectory import TrajectoryWriter
import os
from PIL import Image
import cv2  # from opencv-python
from carmm.utils.povray_render import povray_render, atom_sub


def atoms_to_mp4(atoms, mp4_file, povray=True, generic_projection_settings=None, povray_settings=None, frames_per_second=30,
                 image_rescaling=True, atom_subs=None, keep_temp_files=True, **kwargs):
    """
    A function which takes a list of atoms objects, visualises them in povray with your desired settings and output a .mp4 file.

    Parameters:

    atoms: List of Atoms objects
        Atom images you want to make up your video
    mp4_file: String
        Output file base name (i.e. don't include .mp4)
    povray: Boolean
        Controls whether to call povray to image the atoms objects. Set to False if you already have the images with the
        correct names (e.g. if you've previously run the function, kept the files and want to rerun it with different
        video settings).
    generic_projection_settings: Dictionary
        Settings used by PlottingVariables for automatic rendering
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/utils.py PlottingVariables/__init__ for settings options)
    povray_settings: Dictionary
        Settings used by Povray for automatic rendering
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/pov.py POVRAY/__init__ for settings options)
    frames_per_second: Float
        Speed of switching images (Default is 30 fps)
    image_rescaling: Boolean
        Let the function rescale the images to identical dimensions. If rescaled, the image might move around slightly;
        if not rescaled, you may get transparent sections in your video when displaying images with smaller dimensions.
    atom_subs: List of lists of strings
        Pairs of atomic symbols with the first being changed to the second in all images for clearer visualisation.
        An alternate solution to changing atom colours.
        Tip: Find a second atom with a similar atomic radius to the first with a more distinctive colour
    keep_temp_files: Boolean
        If False, will delete the created .png, .pov and .ini files after use. If True, will
        keep these files and produce a .traj file if any atoms were substituted

    Returns:

    A .mp4 file of the atoms objects, visualised in Povray with the desired settings
    """

    steps = len(atoms)
    digits = len(str(steps - 1))
    indices = [f'%0{digits}d' % x for x in range(steps)]
    filenames = [f'{mp4_file}.{index}.png' for index in indices]

    if povray:
        i = 0
        for image in atoms:
            povray_render(image, output=f'{mp4_file}.{indices[i]}', atom_subs=atom_subs,
                          generic_projection_settings=generic_projection_settings,
                          povray_settings=povray_settings)
            i += 1
            print(f'povray_render: {i}/{steps} ({int(100 * i / steps)}%)')

    mean_width = 0
    mean_height = 0
    if image_rescaling:
        j = 0
        for file in filenames:
            im = Image.open(file)
            w, h = im.size
            mean_width += w
            mean_height += h
            j += 1
            print(f'frame size calculation: {j}/{steps} ({int(100 * j / steps)}%)')

        mean_width = int(mean_width / steps)
        mean_height = int(mean_height / steps)
        k = 0
        for file in filenames:
            im = Image.open(file)
            im_resized = im.resize((mean_width, mean_height), Image.LANCZOS)
            im_resized.save(file, 'png', quality=95)
            k += 1
            print(f'frame resizing: {k}/{steps} ({int(100 * k / steps)}%)')

    try:
        # Set frame from the first image
        frame = cv2.imread(filenames[0])
        height, width, layers = frame.shape

        # filename, fourcc, fps, size
        fourcc = cv2.VideoWriter_fourcc(*'mp4v')
        video = cv2.VideoWriter(f'{mp4_file}.mp4', fourcc, frames_per_second, (width, height))

        # Appending images to video
        for image in filenames:
            video.write(cv2.imread(image))

        # Release the video file
        video.release()
        cv2.destroyAllWindows()
        print('Video Complete!')
    except AttributeError as err:
        print(f'{err}: No images')

    if not keep_temp_files:
        for n in indices:
            os.system(f'rm {mp4_file}.{n}.ini')
            os.system(f'rm {mp4_file}.{n}.pov')
            os.system(f'rm {mp4_file}.{n}.png')

    # For testing purposes only
    return mean_width, mean_height, steps, filenames



def atoms_to_gif(atoms, file, automatic=False, generic_projection_settings=None, povray_settings=None, frames_per_second=30,
                pause_time=0.5, atom_subs=None, gif_options=None, keep_temp_files=True, **kwargs):
    """
    A function which takes a list of atoms objects, visualises them in povray with your desired settings and outputs a .gif file.

    Parameters:

    atoms: List of Atoms objects
        Atom images you want to make up your gif
    file: String
        Output file base name (i.e. don't include .gif)
    automatic: Boolean
        If True, automatically renders images using given settings
        If False, opens ASE GUI to allow the user to manually render the images (**FOLLOW THE GIVEN INSTRUCTIONS**)
    generic_projection_settings: Dictionary
        Settings used by PlottingVariables for automatic rendering
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/utils.py PlottingVariables/__init__ for settings options)
    povray_settings: Dictionary
        Settings used by Povray for automatic rendering
        (see https://gitlab.com/ase/ase/-/blob/master/ase/io/pov.py POVRAY/__init__ for settings options)
    frames_per_second: Float
        Speed of switching images (Default is 30 fps)
    pause_time: Float
        Time (in seconds) to pause on the first and last images (Default is 0.5 seconds)
    atom_subs: List of lists of strings
        Pairs of atomic symbols with the first being changed to the second in
        all images for clearer visualisation. An alternate solution to changing atom colours.
        Tip: Find a second atom with a similar atomic radius to the first with a more distinctive colour
    gif_options: Dictionary of strings
        Settings for the Pillow.Image.save() function. For default setting, don't include
        this parameter. Default options are: "save_all=True, optimize=False, loop=0"
        (see https://pillow.readthedocs.io/en/latest/handbook/image-file-formats.html#gif-saving for full list of options)
    keep_temp_files: Boolean
        If False, will delete the created .png, .pov and .ini files after use. If True, will
        keep these files and produce a .traj file if any atoms were substituted

    Returns:

    A .gif file of the .traj file, visualised in Povray with the desired settings
    """

    steps = len(atoms)

    # Generate the list of povray image filenames
    digits = len(str(steps - 1))
    indices = [f'%0{digits}d' % i for i in range(steps)]
    filenames = [f'{file}.{index}.png' for index in indices]

    if automatic:
        j = 0
        for frame in range(steps):
            frame_atoms = atoms[frame]
            povray_render(frame_atoms, output=f'{file}.{indices[frame]}', view=False, atom_subs=atom_subs,
                          generic_projection_settings=generic_projection_settings, povray_settings=povray_settings)
            j += 1
            print(f'povray_render: {j}/{steps} ({int(100 * j / steps)}%)')
    else:
        if atom_subs is not None:
            for frame in range(steps):
                frame_atoms = atoms[frame]
                atoms[frame] = atom_sub(frame_atoms, atom_subs)
            if keep_temp_files:
                writer = TrajectoryWriter(f'{file}_povray.traj', mode='w')
                for frame in range(steps):
                    writer.write(atoms[frame])
        print(f'***Crucial Steps***\n'
              f'1. In ASE GUI, navigate to Tools -> Render Scene\n'
              f'2. Change "Output basename" to {file}\n'
              f'3. Select "Render all frames"\n'
              f'4. Deselect "Show output window"\n'
              f'5. Change any other settings (e.g. Atomic texture set) as desired')
        view(atoms)
        input('***Press Enter to continue once Povray is finished visualising...***\n')

    gifmaker(file, filenames, frames_per_second, pause_time, gif_options, indices, keep_temp_files)

    print("Happy cooking!")

    # For testing purposes
    return file, steps, atoms, filenames


def gifmaker(file, filenames, frames_per_second, pause_time, gif_options, indices, keep_temp_files):

    duration = [(1 / frames_per_second) * 10**3] * len(filenames)  # Duration in milliseconds

    if pause_time is not None:
        duration[0] = pause_time * 10**3
        duration[-1] = pause_time * 10**3

    # Default gif_options
    if gif_options is None:
        gif_options = {}

    if 'save_all' not in gif_options:
        gif_options['save_all'] = True
    if 'optimize' not in gif_options:
        gif_options['optimize'] = False
    if 'loop' not in gif_options:
        gif_options['loop'] = 0

    # Images
    images = []
    for i in range(len(filenames)):
        try:
            im = Image.open(filenames[i])
            images.append(im)
        except FileNotFoundError as err:
            print(f'{err}')
            break

    try:
        images[0].save(f'{file}.gif', append_images=images[1:], duration=duration, save_all=gif_options['save_all'],
                       optimize=gif_options['optimize'], loop=gif_options['loop'])
    except IndexError as err:
        print(f'{err}: No images')

    # Delete the povray image files if requested
    if not keep_temp_files:
        for index in indices:
            os.system(f'rm {file}.{index}.ini')
            os.system(f'rm {file}.{index}.pov')
            os.system(f'rm {file}.{index}.png')

    # For testing purposes
    return filenames, duration, gif_options
