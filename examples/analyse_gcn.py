# This is testing script for gcn calculator.
from carmm.analyse.gcn import generalised_coordination_number,setup_metal_slab_for_general_coordination_number_calculation 

gcn = setup_metal_slab_for_general_coordination_number_calculation(lattice='fcc', facet=(1,1,1), site='ontop')
assert gcn == 7.5

gcn = setup_metal_slab_for_general_coordination_number_calculation(lattice='fcc', facet=(1,1,1), site='bridge')
assert gcn == 12.5

try:
    gcn = setup_metal_slab_for_general_coordination_number_calculation(lattice='hcp', facet=(1,1,1), site='ontop')
except KeyError:
    "Generalised coordination number works only with FCC and BCC lattice"
