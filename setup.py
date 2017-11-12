#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    seqann Sequence Annotation.
#    Copyright (c) 2017 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#


from setuptools import setup

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'biopython',
    'PyMySQL',
    'six',
    'bson',
    'pytz',
    'numpy'
]

test_requirements = [
    'pytz',
    'biopython',
    'PyMySQL',
    'six',
    'bson',
    'unittest',
    'numpy'
]

setup(
    name='seqann',
    version='0.0.2',
    description="Sequence Annotation",
    long_description=readme + '\n\n' + history,
    author="Mike Halagan",
    author_email='mhalagan@nmdp.org',
    url='https://github.com/nmdp-bioinformatics/seqann',
    packages=[
        'seqann',
        'seqann.models'
    ],
    package_dir={'seqann':
                 'seqann'},
    package_data={'seqann': ['data/*.structure',
                             'data/*.csv', 'data/blast/*',
                             'data/allele_lists/*']},
    install_requires=requirements,
    license="LGPL 3.0",
    zip_safe=False,
    keywords='seqann',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    test_suite='tests',
    tests_require=test_requirements
)
