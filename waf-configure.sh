#!/bin/bash
(./waf distclean && PYTHON=~/Devel/ve26-64/bin/python CYTHON=~/Devel/ve26-64/bin/cython ./waf configure)
