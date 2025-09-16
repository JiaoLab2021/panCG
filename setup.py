import codecs
import os
import re
from setuptools import setup, find_packages
from pathlib import Path


def find_version(info_file):
    info_text = Path(info_file).read_text()
    match = re.search(r"^__version__\s*=\s*['\"]([^'\"]+)['\"]", info_text, re.MULTILINE)
    return match.group(1) if match else "0.0.0"


# these things are needed for the README.md show on pypi (if you dont need delete it)
setup_dir = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(setup_dir, "README.md"), encoding="utf-8") as fh:
    long_description = "\n" + fh.read()

# you need to change all these
VERSION = find_version(os.path.join(setup_dir, "panCG", "lib", "info.py"))
DESCRIPTION = 'An integrative pipeline for family-level super-pangenome analysis across coding and noncoding sequences.'
# LONG_DESCRIPTION = 'xxxxxx'

setup(
    name="panCG",
    version=VERSION,
    author="ltan",
    author_email="lei.tan.bio@outlook.com",
    description=DESCRIPTION,
    long_description_content_type="text/markdown",
    long_description=long_description,
    url='https://github.com/rejo27',
    license='MIT',
    packages=find_packages(),
    entry_points={
        'console_scripts': ['panCG=panCG.panCG:main']
    },
    install_requires=[
        # 'numpy>=2.0.2',
        'numpy<2',
        'pandas>=2.2.2',
        'pyBigWig>=0.3.23',
        'pyranges>=0.1.2',
        'ete3>=3.1.3',
        'biopython>=1.84',
        'networkx>=3.2.1',
        'PyYAML>=6.0.2',
        'scipy>=1.13.1',
        'jcvi>=1.4.21'
    ],
    scripts=[
        'panCG/scripts/CNS.anno.R'
    ],
    keywords=['python', 'panCG', 'windows','mac','linux'],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "Programming Language :: Python :: 3",
        "Operating System :: Unix",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
    ]
)

