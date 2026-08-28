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

<!-- ELARA_X_M07_SCIENTIFIC_DOCS -->
## Scientific guide and paper reproducibility

The repository now includes a model-specific scientific guide covering the
principal mathematics, inputs/outputs, calculation workflow, implementation
mapping and validation context:

- [`docs/SCIENTIFIC_MODEL_GUIDE.md`](docs/SCIENTIFIC_MODEL_GUIDE.md)
- [`docs/PAPER_RESULT_PUBLICATION.md`](docs/PAPER_RESULT_PUBLICATION.md)
- [`results/README.md`](results/README.md)

Paper-linked outputs are accepted only through the publication-safe derived
result policy in `RESULT_PUBLICATION_POLICY.json`.

<!-- ELARA_X_SCIENTIFIC_DOCUMENTATION_START -->
## Scientific documentation

The Elara X repository-enhancement programme provides a consistent scientific guide for this accepted model implementation:

- [Documentation index](docs/README.md)
- [Governing mathematics](docs/MODEL_MATHEMATICS.md)
- [Physical model and process](docs/PHYSICAL_MODEL_AND_PROCESS.md)
- [Inputs, outputs and computation process](docs/INPUTS_OUTPUTS_AND_PROCESS.md)
- [Worked workflow](docs/WORKED_WORKFLOW.md)
- [Validation and reproducibility](docs/VALIDATION_AND_REPRODUCIBILITY.md)
- [Provenance and scientific references](docs/PROVENANCE.md)

Equations in these documents use GitHub-native MathJax syntax.

This enhancement changes documentation only. The accepted scientific source, model resources, licences and validation identity remain unchanged.
<!-- ELARA_X_SCIENTIFIC_DOCUMENTATION_END -->
