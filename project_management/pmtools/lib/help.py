"""help formatting"""
from argparse import HelpFormatter

## I personally don't like the default help formatting output from
## argparse which I think is difficult to read
## Currently identical to RawDescriptionHelpFormatter
class PmHelpFormatter(HelpFormatter):
    def _fill_text(self, text, width, indent):
        return ' '.join([indent + line for line in text.splitlines(True)])

