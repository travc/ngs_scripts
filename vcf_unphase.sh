#!/usr/bin/env bash

awk -F $'\t' '\
BEGIN {OFS = FS}
/^[#]/ {print; next}
{
  for (i = 10; i<=NF; i++) {
    gsub("\\|","/",$i)
  }
print
}
'

# -F $'\t' BEGIN {OFS = FS} # split on tabs
# /^[#]/ {print; next} # just output lines starting with '#' (header)
# the rest replaces '|' with '/' in fields 10+ and prints the line
