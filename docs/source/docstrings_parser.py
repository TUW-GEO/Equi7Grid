# From https://github.com/sphinx-extensions2/sphinx-autodoc2/issues/33
from docutils import nodes
from myst_parser.parsers.sphinx_ import MystParser
from sphinx.ext.napoleon import docstring


class NapoleonParser(MystParser):
    def parse(self, input_string: str, document: nodes.document) -> None:
        parsed_content = str(
            docstring.GoogleDocstring(str(docstring.NumpyDocstring(input_string)))
        )
        return super().parse(parsed_content, document)


Parser = NapoleonParser
