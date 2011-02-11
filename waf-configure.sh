#!/bin/bash
# (./waf distclean && PYTHON=~/Devel/ve26-64/bin/python CYTHON=~/Devel/ve26-64/bin/cython ./waf configure)

(./waf distclean && \
PYTHON=/Library/Frameworks/EPD64.framework/Versions/Current/bin/python \
CYTHON=~/.local/bin/cython \
./waf configure)
