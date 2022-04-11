from rdkit import Chem

class Highlight:
    def __init__(self, indices, text=None, color=None, fill_ring=False, same_color=False):
        """Creates a highlight from a list of atom indices

        Parameters
        ----------
        indices : list[int] or list[list[int]]
            Indices of atoms to highlight
        text : str or None
            Text to match on the molecule's label
        color : str, tuple[float] or None
            Color of the highlight, as a hex code or an RGB tuple. Use `None`
            for automatic color assignment
        fill_ring : bool
            Fill the entire ring or simply highlight the ring atoms
        same_color : bool
            If a nested list is given, assign the same color to all matches
        """
        self.indices = indices
        self.text = text
        self.color = color
        self.fill_ring = fill_ring
        self.same_color = same_color

    def __repr__(self):
        ix = self.indices
        text = self.text
        color = self.color
        fill_ring = self.fill_ring
        return f'<Highlight({ix=}, {text=}, {color=}, {fill_ring=})>'

    @classmethod
    def from_smarts(cls, mol, smarts, text=None, color=None, fill_ring=False,
                    same_color=False, single_match=True):
        """Creates a highlight from a SMARTS string
        
        Parameters
        ----------
        mol : RDKit.Chem.rdchem.Mol
            RDKit mol on which the substructure query is performed
        smarts : str
            SMARTS string for the substructure query
        single_match : bool
            Restrict query matches to a single one
        """
        qmol = Chem.MolFromSmarts(smarts)
        if single_match:
            indices = mol.GetSubstructMatch(qmol)
        else:
            indices = mol.GetSubstructMatches(qmol)
        if not indices:
            raise ValueError(f"No match found for {smarts!r}")
        return cls(indices=indices, text=text, color=color, fill_ring=fill_ring,
                   same_color=same_color)

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, value):
        """RGB tuple (range 0-1) or hex code for color"""
        if isinstance(value, str):
            if value.startswith("#"):
                self._color = value
            else:
                raise ValueError("Please use a hex color code starting with `#`")
        elif isinstance(value, (tuple, list)) and len(value) == 3:
            r, g, b = [int(x*255) for x in value]
            hex = f"#{r:02x}{g:02x}{b:02x}"
            self._color = hex
        elif not value:
            self._color = None
        else:
            raise ValueError("Invalid color format")

    def get_rdkit_color(self):
        """Returns an RGB tuple in the 0-1 range for RDKit"""
        color = self.color.lstrip("#")
        return tuple(int(color[i:i+2], 16)/255 for i in (0, 2, 4))
