#!/bin/bash

(cd test && ./nosetests -v -s test_gsl_interp2d.py test_trace_integrator.py)
