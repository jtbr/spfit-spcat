"""TOML serialization/deserialization for SPFIT/SPCAT input/output structs.

File roles:
  mol.toml         FitInput  — user-authored, read by spfit
  mol.fitted.toml  FitOutput — written by spfit, read by spcat
  mol.dipoles.toml CatInput extras (control + dipoles) — user-authored, read by spcat
  mol.catalog.toml CatOutput — written by spcat

Round-trip workflow:
  mol.toml → spfit mol → mol.fitted.toml
  mol.fitted.toml + mol.dipoles.toml → spcat mol → mol.catalog.toml
"""

import sys
from pathlib import Path
from typing import Union

if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomllib  # type: ignore[no-redef]
    except ImportError:
        try:
            import tomli as tomllib  # type: ignore[no-redef, import]
        except ImportError as exc:
            raise ImportError(
                "Python < 3.11 requires 'tomli' for TOML reading: pip install tomli"
            ) from exc

from ._pickett import (
    FitInput,
    CatInput,
    CalFitOutput,
    CalCatOutput,
    Parameter,
    LineRecord,
    DipoleMoment,
    VibState,
    SpinvOptions,
    DpiOptions,
    EngineOptions,
    EngineKind,
    CatControl,
)

PathLike = Union[str, Path]

# ── Internal helpers ────────────────────────────────────────────────────────


def _vib_to_dict(v: VibState) -> dict:
    d: dict = {}
    if v.knmin != 0:
        d["knmin"] = v.knmin
    if v.knmax != 359:
        d["knmax"] = v.knmax
    if v.stat_weight_axis != 1:
        d["stat_weight_axis"] = v.stat_weight_axis
    if v.iwtpl != 1:
        d["iwtpl"] = v.iwtpl
    if v.iwtmn != 1:
        d["iwtmn"] = v.iwtmn
    if v.vsym != 0.0:
        d["vsym"] = v.vsym
    if v.esym_weight != 99:
        d["esym_weight"] = v.esym_weight
    if v.spin_degeneracies:
        d["spin_degeneracies"] = list(v.spin_degeneracies)
    if v.symmetric_rotor_quanta:
        d["symmetric_rotor_quanta"] = True
    return d


def _vib_from_dict(d: dict) -> VibState:
    v = VibState()
    if "knmin" in d:
        v.knmin = int(d["knmin"])
    if "knmax" in d:
        v.knmax = int(d["knmax"])
    if "stat_weight_axis" in d:
        v.stat_weight_axis = int(d["stat_weight_axis"])
    if "iwtpl" in d:
        v.iwtpl = int(d["iwtpl"])
    if "iwtmn" in d:
        v.iwtmn = int(d["iwtmn"])
    if "vsym" in d:
        v.vsym = float(d["vsym"])
    if "esym_weight" in d:
        v.esym_weight = int(d["esym_weight"])
    if "spin_degeneracies" in d:
        v.spin_degeneracies = [int(x) for x in d["spin_degeneracies"]]
    if "symmetric_rotor_quanta" in d:
        v.symmetric_rotor_quanta = bool(d["symmetric_rotor_quanta"])
    return v


def _spinv_to_dict(s: SpinvOptions) -> dict:
    d: dict = {}
    if s.inclusion_flags != 0:
        d["inclusion_flags"] = s.inclusion_flags
    if s.diag_order != 0:
        d["diag_order"] = s.diag_order
    if s.phase_flags != 0:
        d["phase_flags"] = s.phase_flags
    if s.oblate:
        d["oblate"] = True
    if s.nam_file:
        d["nam_file"] = s.nam_file
    d["vibs"] = [_vib_to_dict(v) for v in s.vibs]
    return d


def _spinv_from_dict(d: dict) -> SpinvOptions:
    s = SpinvOptions()
    if "inclusion_flags" in d:
        s.inclusion_flags = int(d["inclusion_flags"])
    if "diag_order" in d:
        s.diag_order = int(d["diag_order"])
    if "phase_flags" in d:
        s.phase_flags = int(d["phase_flags"])
    if "oblate" in d:
        s.oblate = bool(d["oblate"])
    if "nam_file" in d:
        s.nam_file = str(d["nam_file"])
    if "vibs" in d:
        s.vibs = [_vib_from_dict(v) for v in d["vibs"]]
    return s


def _dpi_to_dict(dpi: DpiOptions) -> dict:
    d: dict = {}
    if dpi.isdgn != 1:
        d["isdgn"] = dpi.isdgn
    if dpi.nvib != 1:
        d["nvib"] = dpi.nvib
    return d


def _dpi_from_dict(d: dict) -> DpiOptions:
    dpi = DpiOptions()
    if "isdgn" in d:
        dpi.isdgn = int(d["isdgn"])
    if "nvib" in d:
        dpi.nvib = int(d["nvib"])
    return dpi


def _engine_options_to_dict(e: EngineOptions) -> dict:
    if e.kind == EngineKind.Dpi:
        d: dict = {"kind": "dpi"}
        sub = _dpi_to_dict(e.dpi)
        if sub:
            d["dpi"] = sub
    else:
        d = {"kind": "spinv"}
        sub = _spinv_to_dict(e.spinv)
        if sub:
            d["spinv"] = sub
    return d


def _engine_options_from_dict(d: dict) -> EngineOptions:
    e = EngineOptions()
    kind_str = str(d.get("kind", "spinv")).lower()
    if kind_str == "dpi":
        e.kind = EngineKind.Dpi
        if "dpi" in d:
            e.dpi = _dpi_from_dict(d["dpi"])
    else:
        e.kind = EngineKind.Spinv
        if "spinv" in d:
            e.spinv = _spinv_from_dict(d["spinv"])
    return e


def _parameter_to_dict(p: Parameter) -> dict:
    d: dict = {"id": p.id, "value": p.value}
    if p.a_priori_error != 1e37:
        d["a_priori_error"] = p.a_priori_error
    if p.fixed:
        d["fixed"] = True
    if p.label:
        d["label"] = p.label
    return d


def _parameter_from_dict(d: dict) -> Parameter:
    p = Parameter()
    p.id = int(d["id"])
    p.value = float(d["value"])
    if "a_priori_error" in d:
        p.a_priori_error = float(d["a_priori_error"])
    if "fixed" in d:
        p.fixed = bool(d["fixed"])
    if "label" in d:
        p.label = str(d["label"])
    return p


def _var_parameter_from_dict(d: dict) -> Parameter:
    """Parse a parameter from mol.fitted.toml; 'error' maps to a_priori_error."""
    p = Parameter()
    p.id = int(d["id"])
    p.value = float(d["value"])
    # 'error' is the fitted 1-sigma; fall back to 'a_priori_error' if present
    if "error" in d:
        p.a_priori_error = float(d["error"])
    elif "a_priori_error" in d:
        p.a_priori_error = float(d["a_priori_error"])
    if "fixed" in d:
        p.fixed = bool(d["fixed"])
    if "label" in d:
        p.label = str(d["label"])
    return p


def _line_record_to_dict(lr: LineRecord) -> dict:
    n = lr.nqn
    qn_used = list(lr.qn[: 2 * n]) if n > 0 else []
    d: dict = {"qn": qn_used, "nqn": n, "freq": lr.freq, "err": lr.err}
    if lr.weight != 1.0:
        d["weight"] = lr.weight
    if lr.blend_tag:
        d["blend_tag"] = lr.blend_tag
    return d


def _line_record_from_dict(d: dict) -> LineRecord:
    lr = LineRecord()
    qn_list = [int(x) for x in d.get("qn", [])]
    lr.qn = qn_list
    lr.nqn = int(d.get("nqn", len(qn_list) // 2))
    lr.freq = float(d["freq"])
    lr.err = float(d.get("err", 1e-7))
    lr.weight = float(d.get("weight", 1.0))
    lr.blend_tag = str(d.get("blend_tag", ""))
    return lr


def _cat_control_to_dict(c: CatControl) -> dict:
    defaults = {
        "iflg": 0, "itag": 999, "qrot": 1000.0, "inblk": 0, "lblk": 0,
        "thrsh": -100.0, "thrsh1": -100.0, "fqmax": 9999.99, "tmq": 300.0, "maxv": -1,
    }
    return {k: getattr(c, k) for k, v in defaults.items() if getattr(c, k) != v}


def _cat_control_from_dict(d: dict) -> CatControl:
    c = CatControl()
    if "iflg"   in d: c.iflg   = int(d["iflg"])
    if "itag"   in d: c.itag   = int(d["itag"])
    if "qrot"   in d: c.qrot   = float(d["qrot"])
    if "inblk"  in d: c.inblk  = int(d["inblk"])
    if "lblk"   in d: c.lblk   = int(d["lblk"])
    if "thrsh"  in d: c.thrsh  = float(d["thrsh"])
    if "thrsh1" in d: c.thrsh1 = float(d["thrsh1"])
    if "fqmax"  in d: c.fqmax  = float(d["fqmax"])
    if "tmq"    in d: c.tmq    = float(d["tmq"])
    if "maxv"   in d: c.maxv   = int(d["maxv"])
    return c


def _dipole_to_dict(dm: DipoleMoment) -> dict:
    d: dict = {"id": dm.id, "value": dm.value}
    if dm.starts_new_component:
        d["starts_new_component"] = True
    if dm.label:
        d["label"] = dm.label
    return d


def _dipole_from_dict(d: dict) -> DipoleMoment:
    dm = DipoleMoment()
    dm.id = int(d["id"])
    dm.value = float(d["value"])
    if "starts_new_component" in d:
        dm.starts_new_component = bool(d["starts_new_component"])
    if "label" in d:
        dm.label = str(d["label"])
    return dm


# ── Public conversion functions ─────────────────────────────────────────────


def fit_input_to_dict(fi: FitInput) -> dict:
    """Convert FitInput → plain dict for TOML serialization (mol.toml)."""
    d: dict = {"title": fi.title}
    if fi.n_iterations != 1:
        d["n_iterations"] = fi.n_iterations
    if fi.marquardt_param != 0.0:
        d["marquardt_param"] = fi.marquardt_param
    if fi.max_obs_calc_err != 1e6:
        d["max_obs_calc_err"] = fi.max_obs_calc_err
    if fi.param_err_scale != 1.0:
        d["param_err_scale"] = fi.param_err_scale
    if fi.freq_scale != 1.0:
        d["freq_scale"] = fi.freq_scale
    if fi.max_lines != 32767:
        d["max_lines"] = fi.max_lines
    if fi.extended_qn:
        d["extended_qn"] = fi.extended_qn
    if fi.nxpar != 0:
        d["nxpar"] = fi.nxpar
    d["engine_options"] = _engine_options_to_dict(fi.engine_options)
    d["parameters"] = [_parameter_to_dict(p) for p in fi.parameters]
    if fi.variance:
        d["variance"] = list(fi.variance)
    # raw_lines (from parse_fit_files legacy round-trip) takes precedence over parsed lines
    if fi.raw_lines:
        d["raw_lines"] = list(fi.raw_lines)
    else:
        d["lines"] = [_line_record_to_dict(lr) for lr in fi.lines]
    return d


def fit_input_from_dict(d: dict) -> FitInput:
    """Build FitInput from a TOML-derived dict (mol.toml)."""
    fi = FitInput()
    fi.title            = str(d.get("title", ""))
    fi.n_iterations     = int(d.get("n_iterations", 1))
    fi.marquardt_param  = float(d.get("marquardt_param", 0.0))
    fi.max_obs_calc_err = float(d.get("max_obs_calc_err", 1e6))
    fi.param_err_scale  = float(d.get("param_err_scale", 1.0))
    fi.freq_scale       = float(d.get("freq_scale", 1.0))
    fi.max_lines        = int(d.get("max_lines", 32767))
    fi.extended_qn      = bool(d.get("extended_qn", False))
    fi.nxpar            = int(d.get("nxpar", 0))
    if "engine_options" in d:
        fi.engine_options = _engine_options_from_dict(d["engine_options"])
    fi.parameters = [_parameter_from_dict(p) for p in d.get("parameters", [])]
    fi.variance   = [float(v) for v in d.get("variance", [])]
    # raw_lines (from legacy round-trip) takes precedence; otherwise use parsed lines
    raw = d.get("raw_lines", [])
    if raw:
        fi.raw_lines = list(raw)
    else:
        fi.lines = [_line_record_from_dict(lr) for lr in d.get("lines", [])]
    return fi


def cat_input_to_dict(ci: CatInput) -> dict:
    """Convert CatInput → plain dict for mol.dipoles.toml (control + dipoles only).

    Engine options and parameters are stored in mol.fitted.toml (written by spfit).
    """
    d: dict = {"title": ci.title}
    ctrl_d = _cat_control_to_dict(ci.control)
    if ctrl_d:
        d["control"] = ctrl_d
    if ci.dipoles:
        d["dipoles"] = [_dipole_to_dict(dm) for dm in ci.dipoles]
    return d


def fit_output_to_dict(out: CalFitOutput, fi: FitInput) -> dict:
    """Build mol.fitted.toml dict from CalFitOutput + FitInput metadata.

    Carries engine_options and parameter labels forward from fi so that spcat
    can read this file without any additional inputs.
    """
    d: dict = {
        "title": fi.title,
        "xsqbest": out.xsqbest,
        "itr": out.itr,
    }
    if fi.extended_qn:
        d["extended_qn"] = True
    d["engine_options"] = _engine_options_to_dict(fi.engine_options)

    params_out = []
    last_independent = -1
    for i, (p, raw_val, raw_err) in enumerate(zip(fi.parameters, out.par, out.erpar)):
        is_dep = p.fixed or p.a_priori_error < 0.0
        if is_dep and last_independent >= 0:
            # out.par[i] is a ratio; reconstruct absolute = ratio × fitted master.
            # Mirrors CalFitIO::writeOutput's dependent-param reconstruction.
            val = raw_val * out.par[last_independent]
            err = abs(raw_val) * out.erpar[last_independent]
        else:
            val, err = raw_val, raw_err
        if not is_dep:
            last_independent = i
        pd: dict = {"id": p.id, "value": val, "error": err}
        if p.fixed:
            pd["fixed"] = True
        if p.label:
            pd["label"] = p.label
        params_out.append(pd)
    d["parameters"] = params_out

    if hasattr(out, "variance") and out.variance:
        d["variance"] = list(out.variance)

    return d


def cat_output_to_dict(out: CalCatOutput) -> dict:
    """Convert CalCatOutput → plain dict for mol.catalog.toml."""
    return {
        "nline": out.nline,
        "temp": list(out.temp),
        "qsum": list(out.qsum),
        "cat_lines": list(out.cat_lines),
    }


# ── File I/O ────────────────────────────────────────────────────────────────


def load_fit_input(path: PathLike) -> FitInput:
    """Read mol.toml → FitInput."""
    with open(path, "rb") as f:
        return fit_input_from_dict(tomllib.load(f))


def load_cat_input(var_path: PathLike, int_path: PathLike) -> CatInput:
    """Read mol.fitted.toml + mol.dipoles.toml → CatInput.

    var_path: path to mol.fitted.toml (written by spfit, contains engine_options +
              fitted parameters + variance).
    int_path: path to mol.dipoles.toml (user-authored, contains control + dipoles).
    """
    with open(var_path, "rb") as f:
        var_d = tomllib.load(f)
    with open(int_path, "rb") as f:
        int_d = tomllib.load(f)

    ci = CatInput()
    ci.title     = str(var_d.get("title", ""))
    ci.int_title = str(int_d.get("title", ""))
    ci.extended_qn = bool(var_d.get("extended_qn", False))
    if "engine_options" in var_d:
        ci.engine_options = _engine_options_from_dict(var_d["engine_options"])
    ci.parameters = [_var_parameter_from_dict(p) for p in var_d.get("parameters", [])]
    ci.variance = [float(v) for v in var_d.get("variance", [])]
    ci.control = _cat_control_from_dict(int_d.get("control", {}))
    ci.dipoles = [_dipole_from_dict(dm) for dm in int_d.get("dipoles", [])]
    return ci


def _toml_scalar(v: object) -> str:
    """Render a scalar value as a TOML literal."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v != v:  # NaN
            return "nan"
        return repr(v)
    if isinstance(v, str):
        escaped = v.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n").replace("\r", "\\r")
        return f'"{escaped}"'
    raise TypeError(f"_toml_scalar: unsupported type {type(v)}")


def _toml_inline_table(d: dict) -> str:
    """Render a dict as a TOML inline table {k = v, ...}."""
    parts = []
    for k, v in d.items():
        if isinstance(v, (dict, list)):
            raise TypeError(f"_toml_inline_table: nested containers not supported for key '{k}'")
        parts.append(f"{k} = {_toml_scalar(v)}")
    return "{" + ", ".join(parts) + "}"


def _toml_dumps(d: dict, prefix: str = "") -> str:
    """Minimal TOML serializer for nested dicts with scalar/list/dict values.

    Supports: scalars, lists of scalars, lists of dicts (→ [[table]] arrays),
    and nested dicts (→ [table] sections). No inline-table nesting.
    """
    lines: list[str] = []
    deferred_sections: list[tuple[str, object]] = []

    for key, val in d.items():
        full_key = f"{prefix}.{key}" if prefix else key
        if isinstance(val, dict):
            deferred_sections.append((full_key, val))
        elif isinstance(val, list) and val and isinstance(val[0], dict):
            deferred_sections.append((full_key, val))
        elif isinstance(val, list):
            if val and isinstance(val[0], dict):
                deferred_sections.append((full_key, val))
            else:
                items = ", ".join(_toml_scalar(item) for item in val)
                lines.append(f"{key} = [{items}]")
        else:
            lines.append(f"{key} = {_toml_scalar(val)}")

    result = "\n".join(lines)
    if lines and deferred_sections:
        result += "\n"

    for full_key, val in deferred_sections:
        if isinstance(val, dict):
            result += f"\n[{full_key}]\n"
            result += _toml_dumps(val, full_key)
        else:  # list of dicts
            for item in val:
                result += f"\n[[{full_key}]]\n"
                result += _toml_dumps(item, full_key)

    return result


def _toml_dump(d: dict, f: object) -> None:
    """Write dict as TOML to a binary file object."""
    f.write(_toml_dumps(d).encode())  # type: ignore[union-attr]


def save_fit_output(out: CalFitOutput, fi: FitInput, path: PathLike) -> None:
    """Write CalFitOutput → mol.fitted.toml."""
    with open(path, "wb") as f:
        _toml_dump(fit_output_to_dict(out, fi), f)


def save_cat_output(out: CalCatOutput, path: PathLike) -> None:
    """Write CalCatOutput → mol.catalog.toml."""
    with open(path, "wb") as f:
        _toml_dump(cat_output_to_dict(out), f)
