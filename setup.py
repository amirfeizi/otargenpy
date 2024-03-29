#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'requests',
    'pandas',
    'matplotlib',
    'seaborn',
    'numpy'
    # Optionally include 'python-Levenshtein' if you want to speed up fuzzywuzzy
]

test_requirements = ['pytest>=3', ]

setup(
    author="Amir Feizi",
    author_email='afeizi@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Interact seamlessly with Open Target Genetics' GraphQL endpoint to query and retrieve tidy data tables, facilitating the analysis of genetic data",
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='otargenpy',
    name='otargenpy',
    packages=find_packages(include=['otargenpy', 'otargenpy.*']),
    package_data={'otargenpy': ['docs/*']},
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/amirfeizi/otargenpy/',
    version='0.2.2',
    zip_safe=False,
)
