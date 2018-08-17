import os
from setuptools import setup, find_packages

# versioning
import versioneer

exepath = 'rstoolbox/bin'
exes = ['check_mutants.py', 'minisilent.py', 'plot_fragments_rmsd.py',
        'regplot_rosetta.py', 'rename_decoys.py'
        ]

setup(
    name='rstoolbox',
    version=versioneer.get_version(),

    description='management and analysis of design populations',
    long_description='rstoolbox is a Python library to visualize, '
                     'analyze and select the designs of interest from a design population. '
                     'It exploits the power of pandas to ease on the selection of decoys of '
                     'interest and provide a direct interface to matplotlib and seaborn plotting.',

    # The project's main homepage.
    url='https://github.com/jaumebonet/RosettaSilentToolbox',

    # Author details
    author='Jaume Bonet',
    author_email='jaume.bonet@gmail.com',

    # Choose your license
    license='MIT',

    # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
    classifiers=[
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Environment :: Console',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Operating System :: MacOS',
        'Operating System :: Unix',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],

    project_urls={
        'Documentation': 'http://jaumebonet.cat/RosettaSilentToolbox',
        'Source': 'https://github.com/jaumebonet/RosettaSilentToolbox/',
        'Tracker': 'https://github.com/jaumebonet/RosettaSilentToolbox/issues',
    },

    platforms='UNIX',
    keywords='development',

    install_requires=['pandas>=0.23', 'pyyaml', 'seaborn', 'libconfig', 'six', 'networkx'],

    packages=find_packages(exclude=['docs', 'demo', 'sphinx-docs']),
    include_package_data=True,
    package_data={
        'rstoolbox': ['analysis/matrices/*',
                      'components/square.ttf',
                      'plot/rama_bgdists/*']
    },
    scripts=[os.path.join(exepath, _) for _ in exes],
    cmdclass=versioneer.get_cmdclass(),
)
