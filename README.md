# MolHighlighter

[![Powered by RDKit](https://img.shields.io/badge/Powered%20by-RDKit-3838ff.svg?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQBAMAAADt3eJSAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAFVBMVEXc3NwUFP8UPP9kZP+MjP+0tP////9ZXZotAAAAAXRSTlMAQObYZgAAAAFiS0dEBmFmuH0AAAAHdElNRQfmAwsPGi+MyC9RAAAAQElEQVQI12NgQABGQUEBMENISUkRLKBsbGwEEhIyBgJFsICLC0iIUdnExcUZwnANQWfApKCK4doRBsKtQFgKAQC5Ww1JEHSEkAAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wMy0xMVQxNToyNjo0NyswMDowMDzr2J4AAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDMtMTFUMTU6MjY6NDcrMDA6MDBNtmAiAAAAAElFTkSuQmCC)](https://www.rdkit.org/)

Fancy colored substructure highlights with labels

<p float="left">
    <img alt="Paired highlight of IUPAC name" height="200" src="assets/iupac_highlight.png"/>
    <img alt="Paired highlight of SMILES string" height="200" src="assets/smiles_highlight.png"/>
</p>

## 🐍 Installation

```python
pip install git+https://github.com/cbouy/molhighlighter.git
```

## 📜 Usage

### Simple substructure highlight

```python
from rdkit import Chem
from molhighlighter import MH, Highlight

mol = Chem.MolFromSmiles("C1=CC(=CC=C1N=NC2=CC=C(C=C2)Br)O")

highlights = [
    Highlight.from_smarts(mol, "c1ccccc1Br", "#93e467", fill_ring=True),
    Highlight.from_smarts(mol, "N=N", "#e36262"),
    Highlight.from_smarts(mol, "c1ccccc1O", "#62d4e3", True),
]
MH(mol, highlights).display()
```

### Substructure highlight paired with label highlighting

```python
from rdkit import Chem
from molhighlighter import LMH, LabelledHighlight

mol = Chem.MolFromSmiles("C1=CC(=CC=C1N=NC2=CC=C(C=C2)Br)O")
name = "4-[(4-bromophenyl)diazenyl]phenol"

highlights = [
    Highlight.from_smarts(mol, "c1ccccc1Br", "#93e467", fill_ring=True),
    Highlight.from_smarts(mol, "N=N", "#e36262"),
    Highlight.from_smarts(mol, "c1ccccc1O", "#62d4e3", True),
]
LMH(name, mol, highlights).display()
```

See the [demo notebook](demo.ipynb) for more info

## ⚖ License

Unless otherwise noted, all files in this directory and all subdirectories are distributed under the Apache License, Version 2.0:
```text
    Copyright 2022 Cédric BOUYSSET

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
```