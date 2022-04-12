import warnings
import uuid
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Draw
from .utils import (Substitution, requires_config,
                    sequential_palette, get_auto_palette)
from .highlight import Highlight
try:
    from IPython.display import display_svg, display_html, Javascript
    from ipywidgets import ColorPicker
except ImportError:
    pass


class MolHighlighter:
    def __init__(self, mol, highlights=None, label=None):
        Chem.GetSSSR(mol)
        self.ring_info = mol.GetRingInfo()
        AllChem.Compute2DCoords(mol)
        self.conformer = mol.GetConformer()
        self.mol = mol
        self.highlights = highlights
        self.label = label
        self._div_id = f"mol_canvas_{uuid.uuid1().hex}"
        self._is_configured = False

    def __repr__(self):
        mol = Chem.MolToSmiles(self.mol)
        label = self.label
        highlights = self.highlights
        return f'MolHighlighter({mol=}, {label=}, {highlights=})'
    
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
        if self.label:
            print(self.label)
        display_svg(svg, raw=True)
        print("Click on the square to pick a color for the highlight")
        return ColorPicker(concise=False, value='#e36262')

    def configure(self, size=(-1, -1), black_font=True, fill_rings=None,
        highlight_font=False, style=None, **moldrawoptions):
        """Configure the highlights

        Parameters
        ----------
        size : tuple[int]
            Size of the molecule image. Use `(-1, -1)` for automatic sizing
        black_font : bool
            Use a black font for the molecule's atoms
        fill_rings : bool or None
            Fill highlighted rings. Leave to `None` to automatically fill
            if rings matching the query were detected
        highlight_font : bool
            Color the label font instead of the label background
        style : str or None
            CSS style for the figure. Not used if label is not set
        moldrawoptions
            See rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
        """
        # check for errors
        if not self.highlights:
            raise AttributeError("Please set the `highlights` attribute")
        if self.label and any(h.text is None for h in self.highlights):
            raise ValueError("Cannot use empty highlight text if label is set")

        # extend highlights if indices is nested list
        extended = []
        for h in self.highlights:
            if isinstance(h.indices[0], (list, tuple)) and not h.same_color:
                ext = [Highlight(ix, h.text, h.color, h.fill_ring)
                       for ix in h.indices]
                extended.extend(ext)
            else:
                extended.append(h)
        self.highlights = extended

        # automatic colors if not set
        if any(h.color is None for h in self.highlights):
            num_hl = len(self.highlights)
            if num_hl > 5:
                palette = get_auto_palette(num_hl)
            else:
                palette = sequential_palette
            for h, color in zip(self.highlights, palette):
                h.color = h.color if h.color else color

        # extend highlights if indices is nested list or text is list
        extended = []
        for h in self.highlights:
            if isinstance(h.indices[0], (list, tuple)):
                ext = [Highlight(ix, h.text, h.color, h.fill_ring)
                       for ix in h.indices]
                extended.extend(ext)
            elif isinstance(h.text, (list, tuple)):
                ext = [Highlight(h.indices, txt, h.color, h.fill_ring)
                       for txt in h.text]
                extended.extend(ext)
            else:
                extended.append(h)
        self.highlights = extended

        # MolDrawOptions defaults
        opts = Draw.MolDrawOptions()
        if black_font:
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
        
        # label
        self.highlight_font = highlight_font

        # html template
        if style:
            self._style = style
        else:
            self._style = ""
        self.html_template = """
        <div id="{div_id}" style="text-align: center; display: inline-block; {style}">
            <style>
                .mhl-highlight {{
                    padding: 1px 1px;
                    border-radius: 6px;
                }}
                .mhl-label {{
                    font-size: 12pt;
                    margin: 6px 0;
                }}
            </style>
            <div class="mhl-label">{label}</div>
            <div>{svg}</div>
        </div>"""
        self._is_configured = True

    @requires_config
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

    def _find_text(self, text, start, occupied):
        index = self.label.find(text, start)
        if index in occupied:
            return self._find_text(text, index + len(text), occupied)
        return index

    def _get_span_element(self, text, color):
        if self.highlight_font:
            return f'<span style="color: {color}">{text}</span>'
        else:
            return ( '<span class="mhl-highlight" '
                    f'style="background-color: {color}">{text}</span>')

    @requires_config
    def generate_label(self):
        """Generate an HTML string of the label with highlights"""
        # sort Highlights with longer (more specific) texts first
        highlights = sorted(self.highlights, reverse=True,
                            key=lambda highlight: len(highlight.text))
        substitutions = []
        highlighted_positions = set()
        for highlight in highlights:
            text = highlight.text
            size = len(text)
            # find start index of text in label
            start = self._find_text(text, 0, highlighted_positions)
            if start < 0:
                warnings.warn(f"{highlight.text!r} unmatched or already found in label")
                continue
            end = start + size
            # create substitution string
            sub = self._get_span_element(text, highlight.color)
            substitutions.append(Substitution(sub, start, end))
            highlighted_positions.update(range(start, end))
        # sort substitutions by order of appearance in label
        substitutions.sort(key=lambda x: x.start)
        n = 0
        label = self.label
        # replace texts by substitutions
        for sub, start, end in substitutions:
            start += n
            end += n
            label = label[:start] + sub + label[end:]
            n += len(sub) - (end - start)
        return label

    @requires_config
    def generate_html(self):
        """Generate the HTML for the labelled figure"""
        label = self.generate_label()
        svg = self.generate_mol_svg()
        return self.html_template.format(div_id=self._div_id, style=self._style,
                                         label=label, svg=svg)

    @requires_config
    def display(self):
        """Displays the highlighted figure"""
        if self.label:
            html = self.generate_html()
            return display_html(html, raw=True)
        svg = self.generate_mol_svg()
        return display_svg(svg, raw=True)

    @requires_config
    def save(self, path=None):
        """Saves the highlighted figure
        
        Parameters
        ----------
        path : str, pathlib.Path or None
            If `label` is None, path to the output SVG file. Else, the HTML view will be
            converted to a PNG and a popup window will open to ask where to save the
            image"""
        # HTML to PNG
        if self.label:
            if path:
                warnings.warn("Path provided but not used for saving")
            self.display()
            script = """
            require.config({
                paths: {
                    html2canvas: 'https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min'
                }
            });
            require(['html2canvas'], function (html2canvas) {
                let div = document.getElementById("%s");
                let options = { backgroundColor: null, imageTimeout: 5000};
                html2canvas(div, options).then(function(canvas) {
                    var link = document.createElement('a');
                    link.download = 'mol_highlight.png';
                    link.href = canvas.toDataURL("image/png").replace("image/png", "image/octet-stream");
                    link.click();
                    link.remove();
                });
            });
            """ % self._div_id
            return Javascript(script)
        # SVG
        if not path:
            raise ValueError("Please set the output path for the SVG file")
        svg = self.generate_mol_svg()
        with open(path, "w") as f:
            f.write(svg)

    def _repr_html_(self):
        if self.label:
            return self.generate_html()

    def _repr_svg_(self):
        return self.generate_mol_svg()
