from paver.easy import *
from paver.setuputils import setup, find_packages
import paver.doctools 

setup(
    name="pm",
    version="0.1.0",
    description = "Project management tools.", 
    author="Per Unneberg",
    url="https://github.com/SciLifeLab/scilifelab.git",
    author_email="per.unneberg@scilifelab.se",
    namespace_packages=["pmtools"],
    packages=find_packages(),
    install_requires = [
        "cement >= 1.9.12",
        ],
    test_suite = 'nose.collector',
    scripts = ['scripts/pm.py'],
    package_dir = {'pmtools':'pmtools'},
    package_data = {'pmtools':['templates/tpl/make/*']}
)


