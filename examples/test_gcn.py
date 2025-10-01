# This is testing script for gcn calculator.
from carmm.analyse.GCN import general_coord_number

gcn = general_coord_number(lattice='fcc', facet=(1,1,1), site='ontop')
print('fcc 111 ontop', gcn)
assert gcn == 7.5
gcn = general_coord_number(lattice='fcc', facet=(1,1,1), site='bridge')
print('fcc 111 bridge', gcn)
assert gcn == 7.5 # No reference for this yet
gcn = general_coord_number(lattice='hcp', facet=(1,1,1), site='ontop')
print('hcp 111 ontop', gcn)
assert gcn == 'gcn failed'
