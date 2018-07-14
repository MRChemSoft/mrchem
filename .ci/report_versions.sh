#!/usr/bin/env bash

set -eu -o pipefail

echo "Report versions of whole tool stack"

for tool in cmake make g++ python; do
    echo ""
    echo "Checking version of $tool:"
    $tool --version 2> /dev/null || echo "$tool not available"
done
