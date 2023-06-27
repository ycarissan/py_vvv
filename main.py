import sys
import rdkit.Chem
import rdkit.Chem.rdMolAlign
import geometry.geometry
import graph_theory.detect_cycle
import numpy as np
import pandas as pd
import plotter.plotter

def main(refname, prbname, a0, b0, c0, a1, b1, c1):
    refMol = rdkit.Chem.MolFromXYZFile(refname)
    prbMol = rdkit.Chem.MolFromXYZFile(prbname)
    jmol_colors_df = pd.read_csv("./jmol_colors.csv")
    print(jmol_colors_df)

    rmsd = rdkit.Chem.rdMolAlign.AlignMol(prbMol, refMol, atomMap=[(a0,a1), (b0,b1), (c0,c1)])
    rdkit.Chem.MolToXYZFile(refMol, "reference.xyz")
    rdkit.Chem.MolToXYZFile(prbMol, "probe.xyz")
    print(rmsd)

    p=plotter.plotter.MyPlotter()

    p.add_ref(refname, jmol_colors_df)
    p.add_prb(prbname, jmol_colors_df)

    p.pl.link_views()
    p.pl.show()

if __name__ == "__main__":
    refname = sys.argv[1]
    prbname = sys.argv[2]
    a0 = int(sys.argv[3])
    b0 = int(sys.argv[4])
    c0 = int(sys.argv[5])
    a1 = int(sys.argv[6])
    b1 = int(sys.argv[7])
    c1 = int(sys.argv[8])
    main(refname, prbname, a0, b0, c0, a1, b1, c1)
