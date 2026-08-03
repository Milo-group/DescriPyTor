from pathlib import Path

from setuptools import setup, find_packages

requirements = [
    line.strip()
    for line in Path('requirements.txt').read_text(encoding='utf-8').splitlines()
    if line.strip() and not line.startswith('#')
]

setup(
    name='DescriPyTor',
    version='0.10005',
    packages=find_packages(exclude=['Getting_started_with_examples*', 'work*']),
    description='Chemical-intuition-based molecular feature extraction and modeling.',
    long_description=Path('README.md').read_text(encoding='utf-8'),
    long_description_content_type='text/markdown',
    author='Eden Specktor',
    author_email='edenpsec@post.bgu.ac.il',
    url='https://github.com/Milo-group/DescriPyTor',
    python_requires='>=3.9',
    install_requires=requirements,
)
