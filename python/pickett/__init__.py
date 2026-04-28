"""
pickett — Python interface to the SPFIT/SPCAT spectroscopy software.

High-level entry points::

    import pickett

    # Fit parameters from .par + .lin files
    fit_out = pickett.fit_files("path/to/molecule")
    print(fit_out.xsqbest, fit_out.par)

    # Generate catalog from .var + .int files (no disk output)
    cat_out = pickett.cat_files("path/to/molecule")
    for line in cat_out.cat_lines[:5]:
        print(line)

Low-level path (access the raw C++ session objects)::

    session = pickett.FitSession("mol.par", "mol.lin", engine="spinv")
    fit_out = session.run()

    session = pickett.CatSession("mol.int", "mol.var", engine="spinv")
    cat_out = session.run()
    q300 = dict(zip(cat_out.temp, cat_out.qsum)).get(300.0)
"""

from ._pickett import (
    CalFitOutput,
    CalCatOutput,
    FitSession,
    CatSession,
    fit_from_files,
    cat_from_files,
    # exceptions
    CalError,
    IoError,
    InputError,
    ValidationError,
    NumericError,
)


def fit_files(base_path: str, engine: str = "spinv") -> CalFitOutput:
    """Fit spectroscopic parameters from *base_path*.par and *base_path*.lin.

    Parameters
    ----------
    base_path:
        Common path prefix, e.g. ``"path/to/co_4/co_4"`` (without extension).
    engine:
        Calculation engine: ``"spinv"`` (default) or ``"dpi"``.

    Returns
    -------
    CalFitOutput
        Fitted parameters (``par``, ``erpar``), RMS (``xsqbest``), and iteration
        count (``itr``).
    """
    return fit_from_files(base_path + ".par", base_path + ".lin", engine)


def cat_files(base_path: str, int_path: str | None = None, engine: str = "spinv") -> CalCatOutput:
    """Generate a spectroscopic catalog from SPCAT input files.

    Reads *base_path*.var (from spfit) and *int_path*.int (dipole/intensity
    file; defaults to *base_path* if not given).  Output is captured in memory
    — no .cat/.egy/.str files are written.

    Parameters
    ----------
    base_path:
        Common path prefix for the .var file.
    int_path:
        Path prefix for the .int file.  Defaults to *base_path*.
    engine:
        Calculation engine: ``"spinv"`` (default) or ``"dpi"``.

    Returns
    -------
    CalCatOutput
        ``cat_lines`` (sorted by frequency), ``egy_lines``, ``str_lines``,
        ``nline``, and partition function tables ``temp``/``qsum``.
    """
    if int_path is None:
        int_path = base_path
    return cat_from_files(int_path + ".int", base_path + ".var", engine)
