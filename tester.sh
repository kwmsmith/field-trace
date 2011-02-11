#!/bin/bash

(cd test && ./nosetests -s test_gsl_interp2d.py test_trace_integrator.py)
