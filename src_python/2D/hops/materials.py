"""Optical constants for refl_map.py: a reader for the refractiveindex.info database
(M. N. Polyanskiy, https://refractiveindex.info, CC0 1.0 public domain) and a curated
library of materials that give interesting HOPS/AWE reflectivity maps.

Refractive indices are complex:  n = n' + i k  (k >= 0 is absorption), matching the
convention of the thesis/paper (n_Ag = 0.05 + 2.275i).

Usage
-----
    from hops.materials import refractive_index, MATERIALS
    refractive_index('Ag', 0.6)                  # curated key, wavelength in microns
    refractive_index('main/Au/nk/Olmon-ev.yml', 1.0, db='/path/to/database/data')
    refractive_index('1.5+0.1j', 0.6)            # a literal index is returned as is

The curated YAML files are bundled in hops/rii_data (copied unchanged from the database).
Any other of the ~3,500 pages can be used by pointing --rii-db (or the environment
variable RII_DB) at the 'database/data' folder of a local copy of
https://github.com/polyanskiy/refractiveindex.info-database  (or the zip from
https://refractiveindex.info/download/database/).
"""
import os
import warnings
from functools import lru_cache

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLED = os.path.join(HERE, 'rii_data')


# ----------------------------------------------------------------------------
# refractiveindex.info YAML reader (formulas 1-9, tabulated n / k / nk)
# ----------------------------------------------------------------------------
def _formula(num, c, lam):
    """Dispersion formulas of refractiveindex.info ('Dispersion formulas.pdf'); lam in um."""
    c = np.concatenate([np.asarray(c, float), np.zeros(17)])
    L = np.asarray(lam, float)
    L2 = L ** 2
    if num == 1:                                   # Sellmeier (preferred)
        n2 = 1 + c[0]
        for i in range(1, 17, 2):
            if c[i] != 0:
                n2 = n2 + c[i] * L2 / (L2 - c[i + 1] ** 2)
        return np.sqrt(n2)
    if num == 2:                                   # Sellmeier-2
        n2 = 1 + c[0]
        for i in range(1, 17, 2):
            if c[i] != 0:
                n2 = n2 + c[i] * L2 / (L2 - c[i + 1])
        return np.sqrt(n2)
    if num == 3:                                   # polynomial
        n2 = c[0] + sum(c[i] * L ** c[i + 1] for i in range(1, 17, 2) if c[i] != 0)
        return np.sqrt(n2)
    if num == 4:                                   # RefractiveIndex.INFO
        n2 = c[0]
        n2 = n2 + c[1] * L ** c[2] / (L2 - c[3] ** c[4])
        n2 = n2 + c[5] * L ** c[6] / (L2 - c[7] ** c[8])
        for i in range(9, 17, 2):
            if c[i] != 0:
                n2 = n2 + c[i] * L ** c[i + 1]
        return np.sqrt(n2)
    if num == 5:                                   # Cauchy
        return c[0] + sum(c[i] * L ** c[i + 1] for i in range(1, 11, 2) if c[i] != 0)
    if num == 6:                                   # gases
        return 1 + c[0] + sum(c[i] / (c[i + 1] - L ** -2) for i in range(1, 11, 2) if c[i] != 0)
    if num == 7:                                   # Herzberger
        return (c[0] + c[1] / (L2 - 0.028) + c[2] * (1 / (L2 - 0.028)) ** 2
                + c[3] * L2 + c[4] * L ** 4 + c[5] * L ** 6)
    if num == 8:                                   # Retro
        A = c[0] + c[1] * L2 / (L2 - c[2]) + c[3] * L2
        return np.sqrt((1 + 2 * A) / (1 - A))
    if num == 9:                                   # Exotic
        n2 = c[0] + c[1] / (L2 - c[2]) + c[3] * (L - c[4]) / ((L - c[4]) ** 2 + c[5])
        return np.sqrt(n2)
    raise ValueError(f'unknown formula {num}')


class RIIMaterial:
    """One refractiveindex.info page (YAML file): callable nk(lambda_um) -> complex."""

    def __init__(self, path):
        import yaml
        with open(path, encoding='utf-8') as fh:
            d = yaml.safe_load(fh)
        self.path = path
        self.references = (d.get('REFERENCES') or '').strip()
        self.comments = (d.get('COMMENTS') or '').strip()
        self.n_part = None       # callable
        self.k_part = None
        lo, hi = [], []
        for D in d.get('DATA', []):
            t = D['type']
            if t.startswith('formula'):
                num = int(t.split()[1])
                coef = [float(v) for v in str(D['coefficients']).split()]
                self.n_part = (lambda L, num=num, coef=coef: np.real(_formula(num, coef, L)))
                r = [float(v) for v in str(D.get('wavelength_range', '0 1e9')).split()]
                lo.append(r[0]); hi.append(r[1])
            elif t.startswith('tabulated'):
                arr = np.array([[float(v) for v in line.split()]
                                for line in str(D['data']).strip().splitlines() if line.strip()])
                arr = arr[np.argsort(arr[:, 0])]
                lam = arr[:, 0]
                lo.append(lam[0]); hi.append(lam[-1])
                if t == 'tabulated nk':
                    self.n_part = (lambda L, a=arr: np.interp(L, a[:, 0], a[:, 1]))
                    self.k_part = (lambda L, a=arr: np.interp(L, a[:, 0], a[:, 2]))
                elif t == 'tabulated n':
                    self.n_part = (lambda L, a=arr: np.interp(L, a[:, 0], a[:, 1]))
                elif t == 'tabulated k':
                    self.k_part = (lambda L, a=arr: np.interp(L, a[:, 0], a[:, 1]))
                # 'tabulated n2' (nonlinear index) is ignored
        if self.n_part is None:
            raise ValueError(f'{path}: no refractive-index data')
        self.range = (max(lo), min(hi)) if lo else (0, np.inf)

    def __call__(self, lam_um, warn=True):
        L = np.asarray(lam_um, float)
        lo, hi = self.range
        if warn and (np.any(L < lo * 0.999) or np.any(L > hi * 1.001)):
            warnings.warn(f'{os.path.basename(self.path)}: wavelength {np.min(L):.4g}-{np.max(L):.4g} um '
                          f'outside the data range {lo:.4g}-{hi:.4g} um (clamped to the range)',
                          RuntimeWarning, stacklevel=2)
        Lc = np.clip(L, lo, hi)
        n = self.n_part(Lc)
        k = self.k_part(Lc) if self.k_part is not None else 0.0 * n
        return n + 1j * k


# ----------------------------------------------------------------------------
# curated library
# ----------------------------------------------------------------------------
# key: (data path inside the database 'data' folder, category, suggested grating period in um
#       (so that the six bands q = 1..6, lambda = period/omega, omega in [1, 7], sample the
#        spectral region where the material is interesting), short description)
M = {}


def _add(key, path, cat, period, note):
    M[key] = dict(path=path, category=cat, period=period, note=note)


# -- superstrates (upper layer, n^u) --------------------------------------------------
_add('air', 'other/mixed gases/air/nk/Ciddor.yml', 'superstrate', 1.0, 'dry air, n = 1.0003')
_add('water', 'main/H2O/nk/Hale.yml', 'superstrate', 1.0, 'water 25 C (0.2-200 um); biosensing, immersion')
_add('ethanol', 'organic/C2H6O - ethanol/nk/Rheims.yml', 'superstrate', 1.0, 'ethanol, n ~ 1.36')
_add('PDMS', 'organic/(C2H6OSi)n - polydimethylsiloxane/nk/Gupta.yml', 'superstrate', 1.0,
     'PDMS elastomer, n ~ 1.41 (stretchable gratings)')
_add('PMMA', 'organic/(C5H8O2)n - poly(methyl methacrylate)/nk/Szczurowski.yml', 'superstrate', 1.0,
     'PMMA resist/cladding, n ~ 1.49')
_add('fused_silica', 'main/SiO2/nk/Malitson.yml', 'superstrate', 1.0, 'fused silica, n ~ 1.46')
_add('BK7', 'specs/schott/optical/N-BK7.yml', 'superstrate', 1.0, 'N-BK7 crown glass, n ~ 1.52 (prism/substrate)')
_add('CaF2', 'main/CaF2/nk/Malitson.yml', 'superstrate', 5.0, 'CaF2, n ~ 1.43, transparent 0.2-9 um')
_add('MgF2', 'main/MgF2/nk/Dodge-o.yml', 'superstrate', 1.0, 'MgF2 (ordinary), n ~ 1.38, low-index coating')
_add('ZnSe', 'main/ZnSe/nk/Connolly.yml', 'superstrate', 10.0, 'ZnSe IR window/ATR prism, n ~ 2.4')
# -- plasmonic (Drude-like, Re eps << 0, small losses) -----------------------------------
_add('Ag', 'main/Ag/nk/Johnson.yml', 'plasmonic metal', 1.0, 'silver, Johnson & Christy 1972 (thesis Fig. 19a)')
_add('Au', 'main/Au/nk/Johnson.yml', 'plasmonic metal', 1.0, 'gold, Johnson & Christy 1972 (thesis Fig. 19b)')
_add('Cu', 'main/Cu/nk/Johnson.yml', 'plasmonic metal', 1.0, 'copper, Johnson & Christy 1972 (thesis Fig. 31a)')
_add('Al', 'main/Al/nk/Rakic.yml', 'plasmonic metal', 0.5, 'aluminium, Rakic 1995: UV plasmonics, interband dip at 0.8 um')
_add('Na', 'main/Na/nk/Smith.yml', 'plasmonic metal', 1.0, 'sodium, Smith 1969: lowest-loss plasmonic metal')
_add('K', 'main/K/nk/Smith.yml', 'plasmonic metal', 1.5, 'potassium, Smith 1969: nearly free-electron')
_add('Mg', 'main/Mg/nk/Palm.yml', 'plasmonic metal', 0.5, 'magnesium, Palm 2018: UV plasmonics')
_add('AuAg50', 'other/alloys/Au-Ag/nk/Rioux-Au50Ag50.yml', 'plasmonic metal', 1.0,
     'Au50Ag50 alloy, Rioux 2014: tunable plasmon between Au and Ag')
_add('Au_IR', 'main/Au/nk/Olmon-ev.yml', 'plasmonic metal', 10.0, 'gold, Olmon 2012 (0.3-25 um): IR Drude metal')
_add('Ag_IR', 'main/Ag/nk/Yang.yml', 'plasmonic metal', 10.0, 'silver, Yang 2015 (0.27-25 um)')
# -- alternative plasmonics: nitrides, TCOs (epsilon-near-zero) ---------------------------
_add('TiN', 'main/TiN/nk/Pfluger.yml', 'alternative plasmonic', 1.0, 'titanium nitride, refractory plasmonic (ENZ ~ 0.5 um)')
_add('ITO', 'other/mixed crystals/In2O3-SnO2/nk/Minenkov-glass.yml', 'ENZ / TCO', 2.0,
     'indium tin oxide on glass, Minenkov 2024: transparent in VIS, metallic in NIR (ENZ ~ 1.3-1.6 um)')
_add('AZO', 'other/doped crystals/Al-ZnO/nk/Shkondin.yml', 'ENZ / TCO', 10.0,
     'Al:ZnO, Shkondin 2017 (2-20 um): ENZ/plasmonic in mid-IR')
# -- lossy / transition metals (large Im n) ---------------------------------------------
_add('W', 'main/W/nk/Werner.yml', 'lossy metal', 1.0, 'tungsten, Werner 2009 (thesis Fig. 20a used Ordal 3.83+2.90i)')
_add('Fe', 'main/Fe/nk/Johnson.yml', 'lossy metal', 1.0, 'iron, Johnson & Christy 1974 (thesis Fig. 20b)')
_add('Co', 'main/Co/nk/Johnson.yml', 'lossy metal', 1.0, 'cobalt, Johnson & Christy 1974 (thesis Fig. 31b)')
_add('Cr', 'main/Cr/nk/Johnson.yml', 'lossy metal', 1.0, 'chromium: strong absorber, adhesion layer')
_add('Ni', 'main/Ni/nk/Johnson.yml', 'lossy metal', 1.0, 'nickel')
_add('Ti', 'main/Ti/nk/Johnson.yml', 'lossy metal', 1.0, 'titanium: n ~ k, near-perfect absorber gratings')
_add('Pt', 'main/Pt/nk/Werner.yml', 'lossy metal', 1.0, 'platinum')
_add('Pd', 'main/Pd/nk/Johnson.yml', 'lossy metal', 1.0, 'palladium (hydrogen sensing)')
_add('Mo', 'main/Mo/nk/Werner.yml', 'lossy metal', 1.0, 'molybdenum')
_add('brass', 'other/alloys/Cu-Zn/nk/Querry-Cu70Zn30.yml', 'lossy metal', 1.0, 'brass Cu70Zn30, Querry 1985')
_add('steel', 'other/alloys/stainless steel/nk/Karlsson-austenitic.yml', 'lossy metal', 1.0,
     'austenitic stainless steel, Karlsson 1982')
# -- phase-change / switchable ------------------------------------------------------------
_add('VO2_cold', 'main/VO2/nk/Beaini-25C.yml', 'phase change', 5.0, 'VO2 at 25 C: insulating phase (0.5-25 um)')
_add('VO2_hot', 'main/VO2/nk/Beaini-100C.yml', 'phase change', 5.0, 'VO2 at 100 C: metallic phase (0.5-25 um)')
# -- semiconductors (high index; absorbing above the band gap) ----------------------------
_add('Si', 'main/Si/nk/Schinke.yml', 'semiconductor', 1.0, 'crystalline silicon, Schinke 2015 (0.25-1.45 um), n ~ 3.5-6.9')
_add('Si_IR', 'main/Si/nk/Franta-25C.yml', 'semiconductor', 5.0, 'silicon, Franta 2017 (0.03-310 um): lossless n = 3.42 in IR')
_add('Ge', 'main/Ge/nk/Nunley.yml', 'semiconductor', 1.0, 'germanium, Nunley 2016 (0.19-2.48 um)')
_add('Ge_IR', 'main/Ge/nk/Amotchkina.yml', 'semiconductor', 10.0, 'germanium, Amotchkina 2020 (0.4-11 um): n = 4.0 in mid-IR')
_add('GaAs', 'main/GaAs/nk/Adachi.yml', 'semiconductor', 1.0, 'gallium arsenide, Adachi 1989 (0.21-12.4 um)')
_add('InP', 'main/InP/nk/Aspnes.yml', 'semiconductor', 0.5, 'indium phosphide, Aspnes 1983')
_add('GaP', 'main/GaP/nk/Aspnes.yml', 'semiconductor', 0.5, 'gallium phosphide, Aspnes 1983: high index, transparent > 0.55 um')
_add('InSb', 'main/InSb/nk/Adachi.yml', 'semiconductor', 5.0, 'indium antimonide, Adachi 1989')
_add('CdTe', 'main/CdTe/nk/Treharne.yml', 'semiconductor', 1.0, 'cadmium telluride (PV absorber)')
_add('MAPbI3', 'other/perovskite/CH3NH3PbI3/nk/Phillips.yml', 'semiconductor', 1.0,
     'methylammonium lead iodide perovskite (PV), Phillips 2015')
_add('MoS2', 'main/MoS2/nk/Song-bulk.yml', 'semiconductor', 1.0, 'bulk MoS2, Song 2019: n ~ 4-5 with exciton resonances')
# -- dielectrics (transparent, index 1.4 - 4) -----------------------------------------
_add('SiO2', 'main/SiO2/nk/Malitson.yml', 'dielectric', 1.0, 'fused silica, Malitson 1965, n ~ 1.46')
_add('sapphire', 'main/Al2O3/nk/Malitson-o.yml', 'dielectric', 1.0, 'sapphire (ordinary), n ~ 1.77')
_add('Si3N4', 'main/Si3N4/nk/Luke.yml', 'dielectric', 1.0, 'silicon nitride, Luke 2015, n ~ 2.0 (photonics)')
_add('TiO2', 'main/TiO2/nk/Devore-o.yml', 'dielectric', 1.0, 'rutile TiO2 (ordinary), n ~ 2.6 (metasurfaces)')
_add('ZnO', 'main/ZnO/nk/Bond-o.yml', 'dielectric', 1.0, 'zinc oxide (ordinary), n ~ 2.0 (thesis Fig. 22 used 2.1054)')
_add('ZnGeP2', 'main/ZnGeP2/nk/Boyd-20C-o.yml', 'dielectric', 5.0,
     'zinc germanium phosphide (ordinary), n ~ 3.2 (thesis Fig. 21 used 3.1874)')
_add('diamond', 'main/C/nk/Phillip.yml', 'dielectric', 1.0, 'diamond, Phillip & Taft 1964, n ~ 2.4')
_add('HfO2', 'main/HfO2/nk/Al-Kuhaili.yml', 'dielectric', 0.5, 'hafnia, n ~ 1.9-2.1 (UV coatings)')
_add('Ta2O5', 'main/Ta2O5/nk/Gao.yml', 'dielectric', 1.0, 'tantala, n ~ 2.1')
_add('LiNbO3', 'main/LiNbO3/nk/Zelmon-o.yml', 'dielectric', 1.0, 'lithium niobate (ordinary), n ~ 2.3 (electro-optic)')
_add('GaN', 'main/GaN/nk/Barker-o.yml', 'dielectric', 1.0, 'gallium nitride (ordinary), n ~ 2.4')
_add('AlN', 'main/AlN/nk/Pastrnak-o.yml', 'dielectric', 0.5, 'aluminium nitride (ordinary), n ~ 2.1')
_add('ZnS', 'main/ZnS/nk/Debenham.yml', 'dielectric', 5.0, 'cubic ZnS, n ~ 2.3 (IR)')
# -- polar dielectrics: phonon polaritons (Re eps < 0 in the Reststrahlen band) ------------
_add('SiC', 'main/SiC/nk/Larruquert.yml', 'phonon polariton', 20.0,
     'SiC, Larruquert 2011 (0.006-132 um): Reststrahlen band 10.3-12.6 um, "silver of the mid-IR"')
_add('SiO2_IR', 'main/SiO2/nk/Popova.yml', 'phonon polariton', 20.0, 'fused silica 7-50 um, Popova 1972: Reststrahlen 8-9.3 um, 20 um')
_add('sapphire_IR', 'main/Al2O3/nk/Querry-o.yml', 'phonon polariton', 25.0,
     'sapphire (ordinary) 0.21-55.6 um, Querry 1985: Reststrahlen 11-17 um')


# -- analytic models (not in the database) ------------------------------------------------
class _Model:
    def __init__(self, fn, rng, ref):
        self.fn, self.range, self.references, self.comments, self.path = fn, rng, ref, '', ref

    def __call__(self, lam_um, warn=True):
        return self.fn(np.asarray(lam_um, float))


def _sic(L):
    # 4H/6H-SiC single Lorentz oscillator (Spitzer 1959; Caldwell et al., Nanophotonics 4, 44 (2015)):
    # eps = eps_inf (1 + (wLO^2 - wTO^2) / (wTO^2 - w^2 - i g w)),  w in cm^-1
    w = 1e4 / L
    eps_inf, wTO, wLO, g = 6.56, 797.0, 970.0, 4.76
    return np.sqrt(eps_inf * (1 + (wLO ** 2 - wTO ** 2) / (wTO ** 2 - w ** 2 - 1j * g * w)))


MODELS = {'SiC': _Model(_sic, (2.0, 50.0), 'Lorentz model, Caldwell et al. 2015: eps_inf=6.56, '
                                            'wTO=797, wLO=970, gamma=4.76 cm^-1')}
M['SiC_film'] = M.pop('SiC')
M['SiC_film']['note'] = 'amorphous SiC thin film, Larruquert 2011 (no sharp Reststrahlen band)'
M['SiC'] = dict(path=None, category='phonon polariton', period=20.0,
                note='crystalline SiC, Lorentz model: Reststrahlen 10.3-12.5 um, "silver of the mid-IR"')

MATERIALS = M


# ----------------------------------------------------------------------------
# public helpers
# ----------------------------------------------------------------------------
def _db_dirs(db=None):
    dirs = [BUNDLED]
    for d in (db, os.environ.get('RII_DB')):
        if d:
            dirs.append(d)
            dirs.append(os.path.join(d, 'data'))
            dirs.append(os.path.join(d, 'database', 'data'))
    return dirs


@lru_cache(maxsize=None)
def load(spec, db=None):
    """Load a curated key or a database path (relative to database/data, or absolute)."""
    if spec in MODELS:
        return MODELS[spec]
    rel = MATERIALS[spec]['path'] if spec in MATERIALS else spec
    if os.path.isabs(rel) and os.path.exists(rel):
        return RIIMaterial(rel)
    for d in _db_dirs(db):
        p = os.path.join(d, rel)
        if os.path.exists(p):
            return RIIMaterial(p)
    raise FileNotFoundError(f"'{spec}' not found (curated keys: see --list-materials; for other pages "
                            f"download the refractiveindex.info database and pass --rii-db)")


def parse_index(spec):
    """Return a complex number if spec is a literal index ('1.1', '0.05+2.275j', '20i'), else None."""
    if isinstance(spec, (int, float, complex)):
        return complex(spec)
    s = str(spec).strip().replace(' ', '').replace('I', 'j').replace('i', 'j')
    try:
        return complex(s)
    except ValueError:
        return None


def refractive_index(spec, lam_um=None, db=None, warn=True):
    """Complex refractive index n + ik of a literal/curated/database material at lam_um (um)."""
    v = parse_index(spec)
    if v is not None:
        return v
    if lam_um is None:
        raise ValueError(f"material '{spec}' needs a wavelength (use --period / --lambda-ref)")
    v = load(spec, db)(lam_um, warn=warn)
    return complex(v) if np.ndim(v) == 0 else np.asarray(v, dtype=complex)


def describe(spec, db=None):
    m = load(spec, db)
    lo, hi = m.range
    info = MATERIALS.get(spec, {})
    return dict(key=spec, range_um=(lo, hi), category=info.get('category', ''), note=info.get('note', ''),
                period=info.get('period'), reference=' '.join(m.references.split())[:160])


def list_materials(lam_um=(0.4, 0.6, 1.0, 1.55, 10.0)):
    """Printable table of the curated library with n at a few wavelengths."""
    rows = []
    for key, info in MATERIALS.items():
        try:
            m = load(key)
        except FileNotFoundError:
            continue
        lo, hi = m.range
        vals = []
        for L in lam_um:
            vals.append(f'{complex(m(L, warn=False)):.3g}' if lo <= L <= hi else '-')
        rows.append((info['category'], key, f'{lo:g}-{hi:g}', info['period'], vals, info['note']))
    return rows
