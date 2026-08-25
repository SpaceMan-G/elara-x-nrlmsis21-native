# NRLMSIS 2.1 — Scientific Model Guide

## Purpose

NRLMSIS 2.1 is an empirical upper-atmosphere specification model. Given time,
geographic position, altitude, solar activity and geomagnetic activity, the
model evaluates neutral temperature and constituent number densities and then
forms the requested aggregate density quantities.

This repository keeps the official coefficient payload external. The public
Python implementation therefore separates **scientific equations and
evaluation logic** from **coefficient/resource authority**.

## Principal inputs

The assembled interface is driven by:

| Quantity | Meaning | Typical unit |
|---|---|---|
| `doy` | day of year | day |
| `utsec` | seconds from 00:00 UT | s |
| `alt` | geodetic altitude | km |
| `glat` | geodetic latitude | degree |
| `glong` | longitude | degree |
| `f107` | daily 10.7 cm solar flux index | sfu |
| `f107a` | centred/mean solar flux index | sfu |
| `ap[0:7]` | geomagnetic activity history | nT-equivalent index |

The official parameter file supplies the empirical coefficients. It is verified
externally before the scientific calculation is permitted to run.

## Mathematical structure

NRLMSIS is not a single closed-form atmosphere equation. It is an empirical
composition of baseline vertical profiles and fitted perturbation functions.

A useful structural representation is

\[
X = X_0 + \sum_j c_j G_j(d,t,\phi,\lambda,F,\bar F,A_p),
\]

where \(X\) may represent a temperature or log-density parameter, \(X_0\) is a
baseline term, \(c_j\) are model coefficients, and \(G_j\) contain the solar,
geomagnetic, seasonal, latitude, diurnal and semidiurnal dependencies.

Above the lower thermosphere, a Bates-like temperature profile is used
conceptually:

\[
T(z) = T_\infty -
\left(T_\infty-T_\ell\right)
\exp[-s\,\xi(z)],
\]

where \(T_\infty\) is exospheric temperature, \(T_\ell\) is a lower-boundary
temperature, \(s\) is a fitted gradient parameter, and \(\xi\) is a
height/geopotential-like coordinate.

For a constituent \(i\), diffusive equilibrium has the form

\[
\frac{d n_i}{dz}
=
-n_i\frac{m_i g}{kT}
-
(1+\alpha_i)\frac{n_i}{T}\frac{dT}{dz},
\]

which integrates to the familiar temperature-ratio and hydrostatic exponential
terms used by the density routines. \(m_i\) is molecular/atomic mass,
\(\alpha_i\) is the thermal-diffusion factor, and \(k\) is Boltzmann's
constant.

The aggregate mass density is formed from constituent number densities:

\[
\rho = \sum_i m_i n_i.
\]

NRLMSIS 2.1 includes the model's accepted nitric-oxide behaviour in the
2.1 state. Exact constituent inclusion follows the public implementation and
legacy-interface semantics rather than an independently reconstructed formula.

## Calculation workflow

```text
verify external parameter authority
        ↓
initialise fitted parameters
        ↓
normalise time / location / solar / geomagnetic inputs
        ↓
evaluate horizontal empirical basis functions
        ↓
evaluate exospheric and local temperature structure
        ↓
evaluate constituent vertical profiles
        ↓
apply lower/upper atmosphere transition logic
        ↓
assemble constituent number densities
        ↓
form aggregate mass density and temperatures
```

## Implementation map

- `parameters.py` — coefficient/state initialisation.
- `horizontal.py` — empirical horizontal, seasonal and activity basis terms.
- `temperature.py` — temperature-profile evaluation.
- `density.py` — constituent density evaluation and vertical integration.
- `model.py` — model-level assembly.
- `legacy_interface.py` — compatibility calculation boundary.
- `api.py` — high-level public initialisation/calculation interface.

## Worked calculation pattern

```python
from elara_x_nrlmsis import initialize, calculate

initialize(resource_file="/path/to/official/resource")

tn, dn = calculate(
    172.0,          # day of year
    43200.0,        # UT seconds
    400.0,          # altitude km
    45.0,           # latitude deg
    10.0,           # longitude deg
    150.0,          # daily F10.7
    155.0,          # mean F10.7
    [4, 4, 4, 4, 4, 4, 4],
)
```

The important reproducibility point is that the returned `tn`/`dn` values are
meaningful only together with the exact repository commit and verified external
parameter-resource identity.

## Validation

The public implementation was accepted only after controlled native/official
equivalence work. The repository safety contract also prohibits redistribution
of official source and official test/resource payloads.

## Reproducible paper-result chain

For paper-linked work, the repository should make the result traceable as:

```text
paper figure/table
        ↓
paper result directory
        ↓
RUN_MANIFEST.json
        ↓
exact model repository commit
        ↓
model inputs and units
        ↓
scientific implementation
        ↓
validation authority / accepted test contract
```

A paper-result directory should therefore record the model name, the exact Git
HEAD used for the calculation, input identities/hashes where redistribution is
permitted, configuration, units, output hashes, and the relationship between
derived plots/tables and the underlying CSV/JSON result.

Only publication-safe derived outputs belong in this repository. Restricted
third-party data, model coefficient payloads that are external by policy,
credentials, private machine paths, private application source, caches, and
temporary files must remain outside the repository.
