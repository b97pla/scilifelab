from datetime import datetime
import codecs

from jinja2 import Environment
from jinja2 import FileSystemLoader


def make_example_note():
    """Make a note with some simple nonsensical data.
    """
    parameters = {
    "formatted_date": "{:%B %d, %Y}".format(datetime.now()),
    "project_name": "A_test",
    "customer_reference": "Some_test",
    "uppnex_project_id": "b2013444",
    "start_date": "000101",
    "FC_id": "SN001_001_AABCD99XX",
    "scilifelab_name": "Test sample",
    "customer_name": "That sample for a test",
    "rounded_read_count": "1",
    "phix_error_rate": "1",
    "avg_quality_score": "1",
    }

    env = Environment(loader=FileSystemLoader(""))
    template = env.get_template("sample_delivery_note.j2")
    with codecs.open("example_sample_delivary_note.html", encoding="utf-8", \
                     mode="w") as note_handle:
        note_handle.write(template.render(parameters))

    print(template.render(parameters))
