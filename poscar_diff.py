import sys

folder_1 = sys.argv[1]
folder_2 = sys.argv[2]

from pymatgen.core import Structure
#from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.structure_matcher import StructureMatcher

s1 = Structure.from_file(f"{folder_1}/POSCAR")
s2 = Structure.from_file(f"{folder_2}/POSCAR")

matcher = StructureMatcher()
if matcher.fit(s1, s2):
    print("same structure")
