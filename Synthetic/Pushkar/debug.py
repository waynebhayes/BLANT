import clean
import os
import sys
script_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, script_path + "/snap3.6/")
sys.path.insert(0, script_path + "/networkx/networkx-2.2/networkx/")
import snap
import networkx

print("sys.path:", sys.path)
print("\n\n")
print("script path:", script_path)
print("\n\n")
print("networkx:", networkx.__dict__)
print("snap:", snap.__dict__)

