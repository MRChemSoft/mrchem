#!/bin/bash
pip install -r doc/requirements.txt
sphinx-build -b html -d doc/doctrees doc doc/html
