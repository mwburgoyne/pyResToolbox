from setuptools import setup, find_packages
#import markdown
import os

ROOT = os.path.abspath(os.path.dirname(__file__))
f = open('C:\\Users\\vinom\OneDrive - Santos\\Work in progress\\ML Learning tidbits\\0-ResEng\\pyResToolbox\\source_dir\\README.md', 'r').read()
long_description= f #markdown.markdown( f.read() )
#long_description = open("README.md").read_text()

setup(
    name = 'pyrestoolbox',
    packages = find_packages(),
    version = '1.1',  # Ideally should be same as your GitHub release tag varsion
    description = 'pyResToolbox - A collection of Reservoir Engineering Utilities',
    long_description= long_description,
    long_description_content_type = 'text/markdown',
    author = 'Mark W. Burgoyne',
    author_email = 'mark.w.burgoyne@gmail.com',
    url = 'https://github.com/mwburgoyne/pyResToolbox',
    #download_url = 'https://github.com/mwburgoyne/pyResToolbox/archive/refs/tags/v1.0.2.tar.gz',
    keywords = ['restoolbox', 'petroleum', 'reservoir'],
    classifiers = [],
    install_requires=[
        'requests',
        'importlib; python_version == "2.6"'
    ],
)