# Elara X NRLMSIS Native

Native Python NRLMSIS 2.1 component developed under the controlled Elara X
atmospheric-model expansion programme.

The translated scientific layers, legacy GTD8D compatibility interface, and
verified external resource resolver are independently validated before the
high-level package assembly is published.

## Parameter resource

The official `msis21.parm` payload is **not bundled** in this repository or in
the public Python package. Runtime initialization requires the externally
supplied official resource with:

- SHA-256: `a322a749f368e73117dd20f3fdcf7389dabc5509f4c27073cc5580999381b508`
- bytes: `536576`

The Stage-9 policy is `EXTERNAL_VERIFIED_RESOURCE_ONLY`.

## High-level interface

The Stage-10 assembled API is designed around:

```python
from elara_x_nrlmsis import initialize, calculate, gtd8d

initialize(resource_file="/path/to/msis21.parm")
tn, dn = calculate(
    172.0, 43200.0, 400.0, 45.0, 10.0,
    150.0, 155.0, [4, 4, 4, 4, 4, 4, 4],
)
```

The assembled calculation boundary fails closed until verified initialization
has completed. Low-level translated modules remain available for controlled
scientific and compatibility use.

## Licence

Use of the NRLMSIS-derived component is governed by
`LICENSE_NRLMSIS21.txt`. The official Fortran source and parameter payload are
not redistributed by this repository.
