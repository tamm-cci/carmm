# This tests the GCN calculator.
from carmm.analyse.gcn import generalised_coordination_number,setup_metal_slab_for_general_coordination_number_calculation 

model, index, cn_max = setup_metal_slab_for_general_coordination_number_calculation(lattice='fcc', facet=(1,1,1), site='ontop')
assert(len(model)) == 1280
assert index  == [476]
assert cn_max == 12

gcn = generalised_coordination_number(model, index, cn_max)
assert gcn == 7.5

model, index, cn_max = setup_metal_slab_for_general_coordination_number_calculation(lattice='fcc', facet=(1,1,1), site='bridge')
assert(len(model)) == 1280
assert index  == [476, 477]
assert cn_max == 12

gcn = generalised_coordination_number(model, index, cn_max)
assert gcn == 12.5

try:
    model, index, cn_max = setup_metal_slab_for_general_coordination_number_calculation(lattice='hcp', facet=(1,1,1), site='ontop')
except KeyError:
    "Model setup only works for FCC and BCC structures"
