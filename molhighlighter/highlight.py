from rdkit import Chem

class Highlight:
    def __init__(self, indices, color, fill_ring=False):
        self.indices = indices
        self.color = color
        self.fill_ring = fill_ring

    @classmethod
    def from_smarts(cls, mol, smarts, color, fill_ring=False):
        qmol = Chem.MolFromSmarts(smarts)
        indices = mol.GetSubstructMatch(qmol)
        if not indices:
            raise ValueError(f"Not match found for {smarts!r}")
        return cls(indices=indices, color=color, fill_ring=fill_ring)

    @property
    def color(self):
        return self._color

    @color.setter
    def color(self, value):
        """RGB tuple or hex code for color"""
        if isinstance(value, str):
            if value.startswith("#"):
                self._color = value
            else:
                raise ValueError("Please use a hex color code starting with `#`")
        elif isinstance(value, (tuple, list)) and len(value) == 3:
            r, g, b = value
            hex = f"#{r:02x}{g:02x}{b:02x}"
            self._color = hex

    def get_rdkit_color(self):
        """Returns an RGB tuple in the 0-1 range for RDKit"""
        color = self.color.lstrip("#")
        return tuple(int(color[i:i+2], 16)/255 for i in (0, 2, 4))


class LabelledHighlight(Highlight):
    def __init__(self, substring, indices, color, fill_ring=False):
        super().__init__(indices=indices, color=color, fill_ring=fill_ring)
        self.substring = substring

    @classmethod
    def from_smarts(cls, mol, smarts, substring, color, fill_ring=False):
        qmol = Chem.MolFromSmarts(smarts)
        indices = mol.GetSubstructMatch(qmol)
        return cls(substring=substring, indices=indices, color=color,
                   fill_ring=fill_ring)
