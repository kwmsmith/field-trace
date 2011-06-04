from field_trace.structure_plonk import structure_plonk

import pylab as pl

def test_structure_plonk():
    arr = structure_plonk(size=512, nstruts=40, amp=2.0, rad_exp=3, core_radius=5, seed=1)
    pl.ion()
    pl.imshow(arr, interpolation='nearest', cmap='hot')
    raw_input('enter to continue')
