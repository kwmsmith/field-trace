
from kaw_analysis import test_vcalc as tv
from kaw_analysis import vcalc
import field_trace

def test_find_nulls():
    test_data = tv.sin_cos_arr(512, 7, 8)
    # test_data *= tv.sin_cos_arr(512, 3, 5)
    d_x = vcalc.cderivative(test_data, 'X_DIR')
    d_y = vcalc.cderivative(test_data, 'Y_DIR')
    nulls = field_trace.find_and_cull_cells(d_x, d_y)
    print len(nulls)
    null_locs = [n.loc for n in nulls]
    X, Y = zip(*null_locs)

if __name__ == '__main__':
    test_find_nulls()
