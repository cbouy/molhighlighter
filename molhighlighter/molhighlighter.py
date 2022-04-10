from collections import namedtuple
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Draw
from IPython.display import display_svg
from ipywidgets import ColorPicker
from .highlight import Highlight


Substitution = namedtuple("Substitution", ["content", "start", "end"])


class MolHighlighter:
    """Highlights substructures and substrings of a molecule and it's
    corresponding label
    """
    highlight_cls = Highlight

    def __init__(self, mol, highlights=None):
        Chem.GetSSSR(mol)
        self.ring_info = mol.GetRingInfo()
        AllChem.Compute2DCoords(mol)
        self.conformer = mol.GetConformer()
        self.mol = mol
        self.highlights = highlights
        self._is_configured = False
    
    def hint(self):
        """Displays the molecule annotated with atom indices and a color picker
        to help setting up the Highlight objects"""
        opts = Draw.MolDrawOptions()
        opts.scalingFactor = 30
        opts.addAtomIndices = True
        opts.annotationFontScale = .7
        d2d = Draw.MolDraw2DSVG(-1, -1)
        d2d.SetDrawOptions(opts)
        d2d.DrawMolecule(self.mol)
        d2d.FinishDrawing()
        svg = d2d.GetDrawingText()
        display_svg(svg, raw=True)
        print("Click on the square to pick a color for the highlight")
        return ColorPicker(concise=False, value='#e36262')

    def configure(self, size=(-1, -1), bw_palette=True, fill_rings=None,
                  **moldrawoptions):
        if not self.highlights:
            raise AttributeError(
                f"Please set the highlights ({self.highlight_cls.__name__}"
                ".highlights) first")
        # MolDrawOptions defaults
        opts = Draw.MolDrawOptions()
        if bw_palette:
            opts.useBWAtomPalette()
        opts.scalingFactor = 28
        opts.highlightBondWidthMultiplier = 16
        opts.highlightRadius = .31
        opts.clearBackground = False
        for prop, value in moldrawoptions.items():
            setattr(opts, prop, value)
        # drawing parameters
        self.size = size
        self.moldrawoptions = opts
        if fill_rings is None:
            self.fill_rings = any(highlight.fill_ring
                                  for highlight in self.highlights)
        else:
            self.fill_rings = fill_rings
        self._is_configured = True

    def generate_mol_svg(self):
        """Creates a drawing of the molecule with highlights"""
        d2d = Draw.MolDraw2DSVG(*self.size)
        d2d.SetDrawOptions(self.moldrawoptions)
        d2d.SetFillPolys(True)
        highlights = {
            "atoms": {},
            "bonds": {}
        }
        mol = self.mol

        # create atom and bond highlight dicts
        for highlight in self.highlights:
            rgb = highlight.get_rdkit_color()
            for i in highlight.indices:
                highlights["atoms"][i] = rgb
                for bond in mol.GetAtomWithIdx(i).GetBonds():
                    bond_atoms = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    if all(atom_ix in highlight.indices for atom_ix in bond_atoms):
                        highlights["bonds"][bond.GetIdx()] = rgb
        
        # draw to set scaling
        d2d.DrawMolecule(
            mol,
            highlightAtoms=list(highlights["atoms"].keys()),
            highlightAtomColors=highlights["atoms"],
            highlightBonds=list(highlights["bonds"].keys()),
            highlightBondColors=highlights["bonds"],
        )
        
        if self.fill_rings:
            # draw filled rings
            for highlight in self.highlights:
                if highlight.fill_ring:
                    rgb = highlight.get_rdkit_color()
                    d2d.SetColour(rgb)
                    ring_atoms = next(ring_atoms
                                      for ring_atoms in self.ring_info.AtomRings()
                                      if all(i in highlight.indices for i in ring_atoms))
                    ps = []
                    for i in ring_atoms:
                        pos = Geometry.Point2D(self.conformer.GetAtomPosition(i))
                        ps.append(pos)
                    d2d.DrawPolygon(ps)

            # draw mol on top of filled rings
            d2d.DrawMolecule(
                mol,
                highlightAtoms=list(highlights["atoms"].keys()),
                highlightAtomColors=highlights["atoms"],
                highlightBonds=list(highlights["bonds"].keys()),
                highlightBondColors=highlights["bonds"],
            )

        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def display(self):
        """Displays the highlighted figure"""
        if not self._is_configured:
            self.configure()
        svg = self.generate_mol_svg()
        return display_svg(svg, raw=True)

    def save(self, path):
        """Saves the highlighted figure
        
        Parameters
        ----------
        path : str or pathlib.Path
            Path to the output SVG file"""
        if not self._is_configured:
            self.configure()
        svg = self.generate_mol_svg()
        with open(path, "w") as f:
            f.write(svg)
