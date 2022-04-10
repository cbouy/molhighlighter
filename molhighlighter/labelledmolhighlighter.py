import warnings
import uuid
from .molhighlighter import Substitution, MolHighlighter
from .highlight import LabelledHighlight
try:
    from IPython.display import display_html, Javascript
except ImportError:
    pass


class LabelledMolHighlighter(MolHighlighter):
    """Highlights substructures of a molecule with substrings of it's
    corresponding label
    """
    highlight_cls = LabelledHighlight

    def __init__(self, label, mol, highlights=None):
        super().__init__(mol=mol, highlights=highlights)
        self.label = label
        self._div_id = f"mol_canvas_{uuid.uuid1().hex}"
        self.style = "text-align: center; display: inline-block;"
        self.html_template = """
        <div id="{div_id}" style="{style}">
            <div style="font-size: 12pt; margin: 10px 0;">{label}</div>
            <div>{svg}</div>
        </div>"""
    
    def hint(self):
        """Displays the molecule annotated with atom indices, its label,
        and a color picker to help setting up the LabelledHighlight objects"""
        print(self.label)
        return super().hint()

    def configure(self, size=(-1, -1), bw_palette=True, fill_rings=None,
                  style="background-color: null",
                  label_style="padding: 4px 1px; border-radius: 6px; background-color:",
                  **moldrawoptions):
        if style:
            self.style += style
        self.label_style = label_style
        super().configure(size=size, bw_palette=bw_palette, fill_rings=fill_rings,
                          **moldrawoptions)

    def _find_substring(self, substring, start, indices):
        index = self.label.find(substring, start)
        if index in indices:
            return self._find_substring(substring, index + len(substring), indices)
        return index

    def generate_label(self):
        """Generate an HTML string of the label with highlights"""
        # sort LabelledHighlights with longer (more specific) substrings first
        highlights = sorted(self.highlights, reverse=True,
                            key=lambda highlight: len(highlight.substring))
        substitutions, starts = [], []
        for highlight in highlights:
            substring = highlight.substring
            size = len(substring)
            # find start index of substring in label
            start = self._find_substring(substring, 0, starts)
            if start < 0:
                warnings.warn(f"No match found in label for {highlight}")
                continue
            end = start + size
            # create substitution string
            sub = f'<span style="{self.label_style}{highlight.color}">{substring}</span>'
            substitutions.append(Substitution(sub, start, end))
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

    def generate_html(self):
        """Generate the HTML for the highlighted figure"""
        if not self._is_configured:
            self.configure()
        label = self.generate_label()
        svg = self.generate_mol_svg()
        return self.html_template.format(div_id=self._div_id, style=self.style,
                                         label=label, svg=svg)

    def display(self):
        """Displays the highlighted figure"""
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
