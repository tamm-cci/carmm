# This is testing script for gcn calculator.
from carmm.analyse.gcn import general_coord_number

gcn = general_coord_number(lattice='fcc', facet=(1,1,1), site='ontop')
assert gcn == 7.5

gcn = general_coord_number(lattice='fcc', facet=(1,1,1), site='bridge')
assert gcn == 12.5 # No reference for this yet

try:
    gcn = general_coord_number(lattice='hcp', facet=(1,1,1), site='ontop')
except KeyError:
    "Generalised coordination number works only with FCC and BCC lattice"
