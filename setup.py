
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'Definition of the Equi7Grid',
    'author': 'Bernhard Bauer-Marschallinger',
    'url': 'https://github.com/bbauerma/Equi7Grid',
    'download_url': 'https://github.com/bbauerma/Equi7Grid',
    'author_email': 'bbm@geo.tuwien.ac.at',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['Equi7Grid'],
    'scripts': [],
    'name': 'Equi7Grid'

}

setup(**config)
