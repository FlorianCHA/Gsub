#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
from pathlib import Path
from setuptools import setup, find_packages
# add for remove error with pip install -e . with pyproject.toml
import site
import sys
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

NAME = 'Gsub'
URL = 'https://github.com/FlorianCHA/Gsub'
CURRENT_PATH = Path(__file__).resolve().parent
VERSION = '2.0.0'

__doc__ = """Gsub is a Graphical user interface that help for GenBank submission.
            For now, this tools take some input file :
                * fasta files which contains all contigs to submit at GenBank
                * template.sbt, a template which create at this link : https://submit.ncbi.nlm.nih.gov/genbank/template/submission/"""

def main():
    setup(
        # Project information
        name=NAME,
        version=VERSION,
        url=URL,
        project_urls={
            'Bug Tracker': f'{URL}/issues',
            'Documentation': URL,
            'Source Code': URL
        },
        download_url=f'{URL}',
        author='''Florian Charriat (CIRAD)''',
        author_email='florian.charriat@cirad.fr',
        description=__doc__.replace("\n",""),
        license='GPLv3',

        # Package information
        packages= ['Gsub'],
        package_data={'Gsub': ['*','image/*','tools/*','tools/tbl2asn.windows/*']},
        include_package_data=True,
        setup_requires=['setuptools_scm'],
        python_requires='>=3.6',
        install_requires=[
            'orffinder',
            'pathlib',
            'Bio',
            'ipython',
            'biopython',
            'colored',
            'gooey==1.0.8.1',
            'argparse',
            'pyinstaller',
            'pyhmmer',
        ],
        # Pypi information
        platforms=['unix', 'linux'],
        keywords=[
            'GenBank submission',
            'python',
            'GUI'
        ],
        classifiers=[
            'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
            'License :: CeCILL-C Free Software License Agreement (CECILL-C)',
            'License :: Free for non-commercial use',
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Developers',
            'Intended Audience :: End Users/Desktop',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Natural Language :: English',
            'Topic :: Scientific/Engineering',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Software Development :: Assemblers'
        ],
        options={
            'bdist_wheel': {'universal': True}
        },

        entry_points={
            f'{NAME}': [f'{NAME} = __init__'],
            'console_scripts': [f'{NAME} = {NAME}.submission_GenBank_UI:IU_parser']},
    )


if __name__ == '__main__':
    main()
