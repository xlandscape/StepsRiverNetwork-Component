project = 'StepsRiverNetwork'
copyright = '2019, knoell Germany GmbH'
author = 'Sebastian Multsch'
version = ''
release = '0.9'
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
pygments_style = 'sphinx'
html_theme = 'classic'



import os
import sys



autoclass_content = 'both'


sys.path.insert(0, os.path.abspath('../..'))
extensions = [
'sphinx.ext.autodoc','sphinx.ext.intersphinx','sphinx.ext.viewcode',]

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
html_static_path = ['_static']
