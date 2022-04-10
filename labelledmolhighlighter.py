import warnings
import uuid
from collections import namedtuple
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Draw
from IPython.display import display_html, display_svg, Javascript
from ipywidgets import ColorPicker


Mapping = namedtuple(
    "Mapping", 
    ["substring", "indices", "color", "fill_ring"],
    defaults=[False])


Subsitution = namedtuple("Substitution", ["content", "start", "end"])


class LabelledMolHighlighter:
    """Highlights substructures and substrings of a molecule and it's
    corresponding label
    """
    def __init__(self, mol, label, mappings=None):
        Chem.GetSSSR(mol)
        self.ring_info = mol.GetRingInfo()
        AllChem.Compute2DCoords(mol)
        self.conformer = mol.GetConformer()
        self.mol = mol
        self.label = label
        self.mappings = mappings
        self._is_configured = False
        self._div_id = f"mol_canvas_{uuid.uuid1().hex}"
        self.html_template = """
        <div id="{div_id}" style="text-align: center; display: inline-block;">
            <div style="font-size: 12pt; margin: 10px 0;">{label}</div>
            <div>{svg}</div>
        </div>"""
    
    def hint(self):
        """Displays the molecule annotated with atom indices, its label,
        and a color picker to help setting up the Mapping objects"""
        opts = Draw.MolDrawOptions()
        opts.scalingFactor = 30
        opts.addAtomIndices = True
        opts.annotationFontScale = .7
        d2d = Draw.MolDraw2DSVG(-1, -1)
        d2d.SetDrawOptions(opts)
        d2d.DrawMolecule(self.mol)
        d2d.FinishDrawing()
        svg = d2d.GetDrawingText()
        print(self.label)
        display_svg(svg, raw=True)
        print("Click on the square to pick a color for the mapping")
        return ColorPicker(concise=False, value='#e36262')

    def configure(self, size=(-1, -1), bw_palette=True, fill_rings=None,
                  label_style="padding: 4px 1px; border-radius: 6px; background-color:",
                  **moldrawoptions):
        if not self.mappings:
            raise AttributeError("Please set the mappings (LabelledMolHighlighter.mappings) first")
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
            self.fill_rings = any(mapping.fill_ring
                                  for mapping in self.mappings)
        else:
            self.fill_rings = fill_rings
        # label parameters
        self.label_style = label_style
        self._is_configured = True

    def hex_to_norm_rgb(self, hexcolor):
        """Convert hex color code to RGB tuple between 0 and 1 for RDKit"""
        color = hexcolor.lstrip("#")
        return tuple(int(color[i:i+2], 16)/255
                     for i in (0, 2, 4))

    def _find_substring(self, substring, start, indices):
        index = self.label.find(substring, start)
        if index in indices:
            return self._find_substring(substring, index + len(substring), indices)
        return index

    def generate_label(self):
        """Generate an HTML string of the label with highlights"""
        # sort Mappings with longer (more specific) substrings first
        mappings = sorted(self.mappings, key=lambda mapping: len(mapping.substring), reverse=True)
        substitutions, starts = [], []
        for mapping in mappings:
            substring = mapping.substring
            size = len(substring)
            # find start index of substring in label
            start = self._find_substring(substring, 0, starts)
            if start < 0:
                warnings.warn(f"No match found in label for {mapping}")
                continue
            end = start + size
            # create substitution string
            sub = f'<span style="{self.label_style}{mapping.color}">{substring}</span>'
            substitutions.append(Subsitution(sub, start, end))
            starts.append(start)
        # sort substitutions by order of appearance in label
        substitutions.sort(key=lambda x: x.start)
        n = 0
        label = self.label
        # replace substrings by substitutions
        for sub, start, end in substitutions:
            start += n
            end += n
            label = label[:start] + sub + label[end:]
            n += len(sub) - (end - start)
        return label

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
        for mapping in self.mappings:
            rgb = self.hex_to_norm_rgb(mapping.color)
            for i in mapping.indices:
                highlights["atoms"][i] = rgb
                for bond in mol.GetAtomWithIdx(i).GetBonds():
                    bond_atoms = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                    if all(atom_ix in mapping.indices for atom_ix in bond_atoms):
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
            for mapping in self.mappings:
                if mapping.fill_ring:
                    rgb = self.hex_to_norm_rgb(mapping.color)
                    d2d.SetColour(rgb)
                    ring_atoms = next(ring_atoms
                                      for ring_atoms in self.ring_info.AtomRings()
                                      if all(i in mapping.indices for i in ring_atoms))
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

    def generate_html(self):
        """Generate the HTML for the highlighted figure"""
        if not self._is_configured:
            self.configure()
        label = self.generate_label()
        svg = self.generate_mol_svg()
        return self.html_template.format(div_id=self._div_id, label=label, svg=svg)

    def display(self):
        """Display the highlighted figure"""
        html = self.generate_html()
        return display_html(html, raw=True)

    def save(self):
        """Displays the highlighted figure and opens a popup to save it"""
        self.display()
        script = """
        require.config({
            paths: {
                html2canvas: 'https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min'
            }
        });
        require(['html2canvas'], function (html2canvas) {
            console.log("doing stuff");
            let div = document.getElementById("%s");
            let options = { backgroundColor: null, imageTimeout: 5000};
            html2canvas(div, options).then(function(canvas) {
                console.log(canvas);
                var link = document.createElement('a');
                link.download = 'mol_highlight.png';
                link.href = canvas.toDataURL("image/png").replace("image/png", "image/octet-stream");
                link.click();
                link.remove();
            });
        });
        """ % self._div_id
        return Javascript(script)
