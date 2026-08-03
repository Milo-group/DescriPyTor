"""Backward-compatible setuptools entry point.

Package configuration lives in ``pyproject.toml``. This file remains for
older tooling that still invokes ``python setup.py`` directly.
"""

from setuptools import setup


if __name__ == "__main__":
    setup()
