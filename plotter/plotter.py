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
        self.actors = []

    def add_mol(self, filename, title, jmol_colors_df, subplot_x=0, subplot_y=0):
        self.pl.subplot(subplot_x, subplot_y)
#        self.pl.enable_mesh_picking(callback=self.oneMeshwasPicked, left_clicking=True)
        self.pl.add_text(title, font_size=30)

        geom = geometry.geometry.Geometry(filename)

        # Adding atoms
        label_points = []
        label_labels = []
        iat=1
        for at in geom.atoms:
            mesh_sphere = pv.Sphere(radius=0.5, center=[at['x'], at['y'], at['z']])
            a=jmol_colors_df[jmol_colors_df.Element==at["label"]]
            try :
                R = a.R.values[0]/255
                G = a.G.values[0]/255
                B = a.B.values[0]/255
            except :
                R=0
                G=255
                B=0
            color=[R, G, B]
            actor = self.pl.add_mesh(mesh_sphere, color=color, show_edges=False)
            label_points.append(mesh_sphere.GetPoint(1))
            label_labels.append(iat)
            iat = iat + 1
            self.actors.append(actor)

        vtk_labels = self.pl.add_point_labels(label_points, label_labels, font_size = 20, always_visible=True)

        # Adding bonds
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
            actor = self.pl.add_mesh(mesh_cylinder, color="white", show_edges=False)

    def oneMeshwasPicked(self, mesh):
        print(mesh)
        print(mesh.getPosition())

