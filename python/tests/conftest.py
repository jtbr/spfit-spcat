"""Change CWD to project root so tests can use spfit_spcat_test_suite/ relative paths."""
import os
from pathlib import Path

# Project root is two levels up from this file (python/tests/conftest.py → repo root)
os.chdir(Path(__file__).parent.parent.parent)
