from distutils.core import setup

setup(
    name = 'pyrestoolbox',
    packages = ['pyrestoolbox'],
    version = 'v1.0.0',  # Ideally should be same as your GitHub release tag varsion
    description = 'pyResToolbox - A collection of Reservoir Engineering Utilities',
    author = 'Mark W. Burgoyne',
    author_email = 'mark.w.burgoyne@gmail.com',
    url = 'https://github.com/mwburgoyne/pyResToolbox',
    download_url = 'download link you saved',
    keywords = ['restoolbox', 'petroleum', 'reservoir'],
    classifiers = [],
    REQUIRED = ['numpy', 'scipy', 'pandas', 'tabulate', 'enum', 'gwr_inversion'],
)
