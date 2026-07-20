#!/bin/sh

cppcheck  --project=compile_commands.json --enable=all --std=c++17 --inline-suppr --suppress=missingIncludeSystem --suppress=checkersReport  -DCPU --check-level=exhaustive  2> 000_cppcheck_advice.txt

# --suppress=unusedFunction
# -v --force
