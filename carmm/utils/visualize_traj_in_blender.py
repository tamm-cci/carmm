default_colors = {
    "H":  (1.0, 1.0, 1.0),
    "C":  (0.05, 0.05, 0.05),
    "N":  (0.1, 0.1, 1.0),
    "O":  (1.0, 0.05, 0.05),
    "F":  (0.56,0.87,0.31),
    "Co": (0.94,0.56,0.62),
    "Ti": (0.74, 0.76, 0.78),
    "Cu": (1.0, 0.85, 0.2),
} # we can add more colors based on Jmol color conventions

vdw_radii = {
    "O":  1.52,
    "F":  1.47,
    "Co": 2.00,
    "Ti": 2.15,
}

def get_color(sym):
    return default_colors.get(sym, (0.6, 0.6, 0.6))

def clear_scene():
    import bpy
    from ase.io import read
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete()

def make_base_sphere(radius, location):
    import bpy
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, location=location)
    base = bpy.context.active_object
    bpy.ops.object.shade_smooth()
    return base

def make_material(base_color):
    import bpy
    mat = bpy.data.materials.new("AtomMat")
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    links = mat.node_tree.links

    # clears any residues and leads to fresh nodes. Nodes important for the color of the material.
    for n in nodes:
        nodes.remove(n)

    out = nodes.new("ShaderNodeOutputMaterial")
    bsdf = nodes.new("ShaderNodeBsdfPrincipled")

    # BSDF, or Bidirectional Scattering Distribution Function,
    # describes how light scatters (reflects and transmits) off a
    # surface from all incoming angles to all outgoing angles, essentially
    # defining a material's appearance.

    # not sure how this affects the scene. Took this from one of the codes in a gitHub project.
    bsdf.location = (0, 0)
    out.location = (300, 0)

    # Set color
    bsdf.inputs["Base Color"].default_value = (*base_color, 1.0)

    # this basically creates a link from the output of the bsdf to input of OutputMaterial
    # I think this is necessary because BSDF tells how the material should behave in
    # presence of light. But, the OutputMaterial helps with actual rendering.
    links.new(bsdf.outputs["BSDF"], out.inputs["Surface"])

    # assigns the same color in viewport and hence no need to go into render mode everytime.
    mat.diffuse_color = (*base_color, 1)

    # in case there are transparent atoms
    mat.blend_method = 'BLEND'
    return mat


def create_structure_in_viewport(atoms, render_cell=False):
    import bpy
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    cell = atoms.cell
    clear_scene() # clearing the scene first
    base_radius = 0.9
    base = make_base_sphere(0.5, (0,0,0))
    atom_objects = []
    for i, (pos, sym, alpha) in enumerate(zip(positions, symbols)):
        radius = vdw_radii[sym]
        obj = base.copy()
        obj.data = base.data.copy()
        bpy.context.collection.objects.link(obj)

        obj.location = pos
        scale = radius / base_radius
        obj.scale = (scale, scale, scale)
        obj.name = f"Atom_{i}_{sym}"

        color = get_color(sym)
        mat = make_material(color)
        obj.data.materials.append(mat)

        atom_objects.append(obj)

    if render_cell:
        make_cell(cell=cell)

    # Hide the base sphere as it is not necessary to view it
    base.hide_viewport = True
    base.hide_render = True
    return atom_objects

# lets create a cell for the structure (if it has one!!)
# the idea is to give the vectors that form the corners of the box and then
# connect the corners using edges (this could be used universally to create a box.

def make_cell(cell):
    import bpy
    from mathutils import Vector
    corners = [
        Vector([0, 0, 0]),
        Vector(cell[0]),
        Vector(cell[1]),
        Vector(cell[2]),
        Vector(cell[0] + cell[1]),
        Vector(cell[0] + cell[2]),
        Vector(cell[1] + cell[2]),
        Vector(cell[0] + cell[1] + cell[2])
    ]

    edges = [
        (0,1),(0,2),(0,3),
        (1,4),(1,5),
        (2,4),(2,6),
        (3,5),(3,6),
        (4,7),(5,7),(6,7)
    ]

    mesh = bpy.data.meshes.new("CellBox") # "CellBox" is not arbitrary and blender looks for this keyword
    mesh.from_pydata(corners, edges, [])
    cell_obj = bpy.data.objects.new("CellBox", mesh)
    bpy.context.collection.objects.link(cell_obj)

    cell_obj.display_type = 'WIRE'
    cell_obj.show_wire = True

    return cell_obj
