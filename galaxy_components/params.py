"""Parameter files.

The native format is **TOML** (standard library, no dependency): one table per group, real
lists, real booleans, comments.  The legacy Fortran namelist format (``&HALO ... /``) is still
read, for the old ``parameters.txt``, if ``f90nml`` is installed.  The format is chosen by
the file suffix: ``.toml`` is TOML, anything else is read as a namelist.

    [halo]
    m200 = 100.0
    c200 = 10.0

    [nbody]
    npart_comp = [20000]          # lists may be shorter than 5; the rest keep their defaults
    form_comp = ["nfw"]

TOML files are validated strictly (an unknown group or field is an error, so a typo cannot be
silently ignored); namelists stay lenient because the Fortran groups carry fields this code
does not use.  ``galaxy-components convert old.txt -o new.toml`` converts a namelist.

The compile-time options of the Fortran code (-DNFW_HALO, -DHERNQUIST_BULGE, -DGalIC, ...)
are runtime options in the ``[model]`` group.  Defaults reproduce the makefile as shipped.
"""
import tomllib
import warnings
from dataclasses import dataclass, field, fields
from pathlib import Path

NCOMP_MAX = 5


@dataclass
class Params:
    # [halo]
    m200: float = 0.0
    v200: float = 0.0
    c200: float = 10.0
    qflat: float = 1.0
    rcore: float = 0.0
    # [disc]
    mdisc: float = 0.0
    sigma0: float = 0.0
    fdisc: float = 1.0
    mstar: float = 0.0
    fhole: float = 0.0
    fdref: float = 0.0
    # [bulge]
    mbulge: float = 0.0
    fbulge: float = 0.0
    # [gas_disc]
    temp: float = 1e4
    ndisc: int = 0
    # [gas_halo]
    fbaryon: float = 0.16
    nhalo: int = 0
    v_radial: float = 0.0
    rtrunc: float = 1.0
    form: str = "nfw"
    # [black_hole]
    mbh: float = 0.0
    # [nbody]
    ncomp: int = 1
    npart_comp: list = field(default_factory=lambda: [0] * NCOMP_MAX)
    fmass_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    rscale_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    raniso_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)  # Osipkov-Merritt radius / rs (0 = off)
    beta_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)    # constant anisotropy beta (0 = isotropic)
    rtrunc_comp: list = field(default_factory=lambda: [0.0] * NCOMP_MAX)
    form_comp: list = field(default_factory=lambda: [""] * NCOMP_MAX)
    istaper: bool = False
    # [output]
    disc_file: str = "gas_disc"
    halo_file: str = "gas_halo"
    nbody_file: str = "nbody"
    snap_format: int = 3
    # [model]  (Python only; replaces the Fortran -D flags)
    halo_type: str = "hernquist"  # 'hernquist' | 'nfw' | 'logarithmic'
    galic: bool = True  # GalIC conventions (total halo mass == mhalo)
    halo_zero_edge: bool = False  # hydrostatic gas out to rtrunc*r200 with P=0 there (Donnert 2013)
    halo_core: float = 0.0  # gas entropy floor, in units of the K&S entropy at r200 (0 = pure K&S)
    # gas-disc mesh (hardcoded in mk_gas_disc.f90)
    nrbin: int = 150
    nzbin: int = 500
    log_zmin: float = -15.0
    log_zmax: float = 1.3
    snap_z: bool = False  # snap vertical positions to the z table, as Fortran


# Which group each field lives in (the namelist groups of the Fortran code, plus ``model``).
GROUP_FIELDS = {
    "halo": ("m200", "v200", "c200", "qflat", "rcore"),
    "disc": ("mdisc", "sigma0", "fdisc", "mstar", "fhole", "fdref"),
    "bulge": ("mbulge", "fbulge"),
    "gas_disc": ("temp", "ndisc"),
    "gas_halo": ("fbaryon", "nhalo", "v_radial", "rtrunc", "form"),
    "black_hole": ("mbh",),
    "nbody": ("ncomp", "npart_comp", "fmass_comp", "rscale_comp", "raniso_comp", "beta_comp",
              "rtrunc_comp", "form_comp", "istaper"),
    "output": ("disc_file", "halo_file", "nbody_file", "snap_format"),
    "model": ("halo_type", "galic", "halo_zero_edge", "halo_core",
              "nrbin", "nzbin", "log_zmin", "log_zmax", "snap_z"),
}
# Fields of the Fortran namelists that this code never used or does not implement.  They are not
# silently ignored: TOML raises, a legacy namelist warns (and the value is dropped).
REMOVED = {
    "sigma200": ("never used by the Fortran or this code", 0.0),
    "lambda_b": ("gas-halo spin was commented out in the Fortran and is not implemented", 0.0),
    "fdecay": ("never used; the exponential taper length is fixed by the slope at rtrunc_comp", 3000.0),
    "glass_file": ("glass sampling is not implemented; the gas is Poisson-sampled", ""),
    "ispoisson": ("only Poisson sampling is implemented", True),
    "tabfile": ("tabulated gas haloes are not implemented", ""),
    "istabulated": ("tabulated gas haloes are not implemented", False),
    "exponential_disc": ("only the exponential disc potential is implemented", True),
}
_LIST_FIELDS = ("npart_comp", "fmass_comp", "rscale_comp", "raniso_comp", "beta_comp", "rtrunc_comp",
                "form_comp")
_TYPES = {f.name: f.type for f in fields(Params)}
_FIELD_GROUP = {name: g for g, names in GROUP_FIELDS.items() for name in names}


def read_params(path) -> Params:
    """Read a parameter file: ``.toml`` as TOML, anything else as a Fortran namelist."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"parameter file does not exist: {path}")
    p = _read_toml(path) if path.suffix.lower() == ".toml" else _read_namelist(path)
    p.halo_type = p.halo_type.lower()
    return p


# ------------------------------------------------------------------ TOML
def _coerce(name, value, where):
    """Check a TOML value against the field's type (ints are accepted for floats)."""
    want = _TYPES[name]
    ok = (isinstance(value, bool) if want is bool else
          isinstance(value, (int, float)) and not isinstance(value, bool) if want is float else
          isinstance(value, int) and not isinstance(value, bool) if want is int else
          isinstance(value, str) if want is str else True)
    if not ok:
        raise ValueError(f"{where}: '{name}' must be {want.__name__}, got {value!r}")
    return float(value) if want is float else value


def _coerce_list(name, values, where):
    if not isinstance(values, list):
        raise ValueError(f"{where}: '{name}' must be a list, got {values!r}")
    if len(values) > NCOMP_MAX:
        raise ValueError(f"{where}: '{name}' has {len(values)} entries; at most {NCOMP_MAX} components")
    cur = list(getattr(Params(), name))
    for i, v in enumerate(values):
        if name == "form_comp":
            if not isinstance(v, str):
                raise ValueError(f"{where}: '{name}' entries must be strings, got {v!r}")
            cur[i] = v.strip().lower()
        elif name == "npart_comp":
            if isinstance(v, bool) or not isinstance(v, int):
                raise ValueError(f"{where}: '{name}' entries must be integers, got {v!r}")
            cur[i] = v
        else:
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise ValueError(f"{where}: '{name}' entries must be numbers, got {v!r}")
            cur[i] = float(v)
    return cur


def _read_toml(path: Path) -> Params:
    with open(path, "rb") as fh:
        try:
            data = tomllib.load(fh)
        except tomllib.TOMLDecodeError as e:
            raise ValueError(f"{path}: invalid TOML: {e}") from None
    p = Params()
    for group, table in data.items():
        for key in (table if isinstance(table, dict) else ()):
            if key in REMOVED:
                raise ValueError(f"{path}: '{key}' in [{group}] is no longer supported: {REMOVED[key][0]}")
        if group not in GROUP_FIELDS:
            raise ValueError(f"{path}: unknown table [{group}]; valid tables: {', '.join(GROUP_FIELDS)}")
        if not isinstance(table, dict):
            raise ValueError(f"{path}: [{group}] must be a table")
        for key, val in table.items():
            if key not in GROUP_FIELDS[group]:
                where = f" (it belongs in [{_FIELD_GROUP[key]}])" if key in _FIELD_GROUP else ""
                raise ValueError(f"{path}: unknown field '{key}' in [{group}]{where}; "
                                 f"valid fields: {', '.join(GROUP_FIELDS[group])}")
            setattr(p, key, _coerce_list(key, val, f"{path} [{group}]") if key in _LIST_FIELDS
                    else _coerce(key, val, f"{path} [{group}]"))
    return p


# ------------------------------------------------------------------ namelist (legacy)
def _read_namelist(path: Path) -> Params:
    try:
        import f90nml
    except ImportError:
        raise ImportError(f"{path}: reading a Fortran namelist needs the optional 'f90nml' package "
                          "(pip install f90nml), or convert it once to TOML with a machine that has it: "
                          "galaxy-components convert FILE -o FILE.toml") from None
    nml = f90nml.read(str(path))
    p = Params()
    for group in nml:
        for key, val in nml[group].items():
            key = key.lower()
            if key in REMOVED:
                if val != REMOVED[key][1] and not (isinstance(val, str) and val.strip() == REMOVED[key][1]):
                    warnings.warn(f"{path}: '{key}' = {val!r} is ignored: {REMOVED[key][0]}", stacklevel=3)
                continue
            if not hasattr(p, key):
                continue
            if key in _LIST_FIELDS:
                cur = list(getattr(p, key))
                vals = val if isinstance(val, list) else [val]
                for i, v in enumerate(vals[:NCOMP_MAX]):
                    if v is not None:
                        cur[i] = v.strip().lower() if isinstance(v, str) else v
                setattr(p, key, cur)
            else:
                setattr(p, key, val.strip() if isinstance(val, str) else val)
    return p


# ------------------------------------------------------------------ TOML writer / converter
def _toml_value(v):
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        return repr(v)
    if isinstance(v, str):
        return '"' + v.replace("\\", "\\\\").replace('"', '\\"') + '"'
    if isinstance(v, list):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    raise TypeError(f"cannot write {v!r} as TOML")


def to_toml(p: Params, header="", all_fields=False) -> str:
    """Serialise ``p`` as TOML.  By default only fields that differ from the defaults are written."""
    default = Params()
    out = [f"# {line}" for line in header.splitlines()] + ([""] if header else [])
    for group, names in GROUP_FIELDS.items():
        lines = []
        for name in names:
            v, d = getattr(p, name), getattr(default, name)
            if not all_fields and v == d:
                continue
            if name in _LIST_FIELDS:      # drop trailing entries that equal the default, keep ncomp
                keep = max(p.ncomp, 1)
                v = list(v)
                while len(v) > keep and v[-1] == d[len(v) - 1]:
                    v.pop()
            lines.append(f"{name} = {_toml_value(v)}")
        if lines:
            out += [f"[{group}]"] + lines + [""]
    return "\n".join(out).rstrip() + "\n"


def output_name(base, snap_format):
    suffix = {1: ".snp1.gdt", 2: ".snp2.gdt", 3: ".gdt.hdf5", 4: ".swft.hdf5"}
    return base + suffix.get(snap_format if snap_format <= 4 else 1, "")
