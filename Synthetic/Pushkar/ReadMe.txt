Notes from WH: This is Pushkar's code that computes a bunch of standard measures for a given network. It probably needs
to be updated to a more recent version of Python... but it would be MUCH better if it were re-written in C, C++, or awk.

--|Notes from AN: This code has since been converted to C++|--
Setup:
Install 7zip if not already on system
Extract SNAP-6.0.7z with 7zip

Usage:
Run with command g++ -ISnap-6.0/glib-core main.cpp
------------------------------------------------------------------

-- Original Readme from Pushkar Gupta, around 2019.--
Runs with Python3.
May need NetworkX library, and SNAP network database.

Tested with python 3.4, 3.5, 3.6
External librares used: SNAP (Python 3.6) - which is an experimental version

## Files
main.py - The actual tool
clean.py - Converts graph with string node-ids to graph with integer node-ids (NO need to call this)

## Usage
python main.py [-e eigenVals] inputFile  # graph can have strings/int node-ids; skip -e param to compute all eigenvalues

