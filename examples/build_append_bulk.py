#!/usr/bin/env python3

'''
Testing appending bulk beneath surfaces. Useful for building
surface QM/MM models
'''

def test_append_bulk():
    from carmm.build.append_bulk import append_bulk
    from ase.build import bulk, surface
    from ase.io import read

    mgo = bulk("MgO", crystalstructure="rocksalt", a=4.212, cubic=True)
    slab=surface(mgo, (1,0,0),3, vacuum=5.0)

    apnd = append_bulk(surface=slab, bulk=bulk, nlayers=3, vac=10.0, output=None)

    apnd.write("atoms_test_append.traj")
    size=apnd.positions.size
    shape=apnd.positions.shape
    assert(size == 144)
    assert(shape == (48,3))

    append_bulk(surface="data/mgo/mgo_slab.traj", bulk="data/mgo/mgo_bulk.traj", nlayers=None, vac=10.0, output="data/mgo/file_test_append.traj")

    test=read("data/mgo/file_test_append.traj")
    size=test.positions.size 
    shape=test.positions.shape
    assert(size == 96)
    assert(shape == (32,3))

test_append_bulk()
