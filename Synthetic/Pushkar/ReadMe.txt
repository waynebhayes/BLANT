Notes from WH: This is Pushkar's code that computes a bunch of standard measures for a given network. It probably needs
to be updated to a more recent version of Python... but it would be MUCH better if it were re-written in C, C++, or awk.

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

