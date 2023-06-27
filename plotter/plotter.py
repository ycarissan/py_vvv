import geometry.geometry
import graph_theory.detect_cycle
import numpy as np

try :
    import pyvista as pv
except ModuleNotFoundError as error:
    pyvista = None
    print("pyvista module not found")

class MyPlotter():
    def __init__(self):
        self.pl = pv.Plotter(shape=(1,2))

    def add_ref(self, refname, jmol_colors_df):
        self.pl.subplot(0, 0)
        self.pl.add_text("Reference", font_size=30)
        spheres, colors, cylinders = get_spheres_colors_and_cylinders(refname, jmol_colors_df)
        for sphere,color in zip(spheres, colors):
            self.pl.add_mesh(sphere, color=color, show_edges=False)
        for cyl in cylinders:
            self.pl.add_mesh(cyl, color="white", show_edges=False)

    def add_prb(self, prbname, jmol_colors_df):
        self.pl.subplot(0, 1)
        self.pl.add_text("Probe", font_size=30)
        spheres, colors, cylinders = get_spheres_colors_and_cylinders(prbname, jmol_colors_df)
        for sphere,color in zip(spheres, colors):
            self.pl.add_mesh(sphere, color=color, show_edges=False)
        for cyl in cylinders:
            self.pl.add_mesh(cyl, color="white", show_edges=False)

def get_spheres_colors_and_cylinders(filename, jmol_colors_df):
    geom = geometry.geometry.Geometry(filename)
    spheres = []
    cylinders = []
    colors = []
    for at in geom.atoms:
        mesh_sphere = pv.Sphere(radius=0.5, center=[at['x'], at['y'], at['z']])
        a=jmol_colors_df[jmol_colors_df.Element==at["label"]]
        try :
            R = a.R.values[0]/255
            G = a.G.values[0]/255
            B = a.B.values[0]/255
        except :
#            print(a)
#            print(at)
            R=0
            G=255
            B=0
        color=[R, G, B]
        spheres.append(mesh_sphere)
        colors.append(color)
    molecularGraph = graph_theory.detect_cycle.MolecularGraph(filename)
    for e in molecularGraph.getEdges():
        idx1 = e.GetBeginAtomIdx()
        idx2 = e.GetEndAtomIdx()
        at1 = geom.getAtom(idx1)
        at2 = geom.getAtom(idx2)
        pos1 = np.asarray([at1['x'], at1['y'], at1['z']])
        pos2 = np.asarray([at2['x'], at2['y'], at2['z']])
        vect_bond = pos2 - pos1
        middle_bond = 0.5 * (pos1 + pos2)

        mesh_cylinder = pv.Cylinder(center=middle_bond, direction=vect_bond, radius=.1, height=np.linalg.norm(vect_bond))
        cylinders.append(mesh_cylinder)
    return spheres, colors, cylinders
