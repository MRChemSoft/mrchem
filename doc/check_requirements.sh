#!/bin/bash
requirements="yaml matplotlib recommonmark sphinx breathe"
for r in $requirements
do
    python -c "import re, $r; print(re.compile('/__init__.py.*').sub('',$r.__file__))"
done
