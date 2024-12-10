#!/bin/sh
cd canon_maps && ls | hawk '/[A-Z]/{printf "ln -sf %s %s\n", $0, tolower($0)}' | sh
