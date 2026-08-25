# Elara X NRLMSIS 2.1 Native

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
`nrlmsis2.1_license.txt`. The official Fortran source and parameter payload are
not redistributed by this repository.

## Repository identity

The M01 independent audit conclusively classified the historical repository as NRLMSIS 2.1 at base commit `1ca8b3cf1a29ad6b0d3b328759ab9381e2e23774` and tree `77c2d58866be0498900f04d9774da3e4657f2ab4`. The eight scientific translation modules remain byte-identical to the accepted Elara X NRLMSIS 2.1 core.

This is a development identity (`0.1.0.dev0`). Final paper-linked tags and releases are deferred until the complete atmospheric-model freeze.
