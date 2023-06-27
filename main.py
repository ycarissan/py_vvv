import sys
import rdkit.Chem
import rdkit.Chem.rdMolAlign

refname = sys.argv[1]
prbname = sys.argv[2]
a0 = int(sys.argv[3])
b0 = int(sys.argv[4])
c0 = int(sys.argv[5])
a1 = int(sys.argv[6])
b1 = int(sys.argv[7])
c1 = int(sys.argv[8])
refMol = rdkit.Chem.MolFromXYZFile(refname)
prbMol = rdkit.Chem.MolFromXYZFile(prbname)
rmsd = rdkit.Chem.rdMolAlign.AlignMol(prbMol, refMol, atomMap=[(a0,a1), (b0,b1), (c0,c1)])
rdkit.Chem.MolToXYZFile(refMol, "reference.xyz")
rdkit.Chem.MolToXYZFile(prbMol, "probe.xyz")
print(rmsd)
