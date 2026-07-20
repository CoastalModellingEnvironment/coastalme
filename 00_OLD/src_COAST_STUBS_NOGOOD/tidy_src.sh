#!/bin/bash

if [ -z "$1" ]; then
  echo "Usage: $0 <filename>"
  exit 1
fi

echo "Formatting $1" && clang-format -i "$1"

