"""Compatibility namespace for code written before the rename to vamp.

New code should import :mod:`vamp`.  Submodules remain discoverable from the
same source path so existing scripts can migrate without losing archive access.
"""
from vamp import __path__, __version__
