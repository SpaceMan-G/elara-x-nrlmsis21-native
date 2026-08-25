# Paper-result publication and reproducibility

## Purpose

This repository is designed to hold **publication-safe model outputs** produced
for scientific papers, reports and validation studies.

Results are stored under:

```text
results/paper/<paper-id>/<execution-id>/
```

Each execution directory should contain at minimum:

- `RUN_MANIFEST.json`;
- `CHECKSUMS.sha256`;
- one or more derived result artefacts such as CSV/JSON tables or publication
  figures.

## What may be published

Suitable examples include:

- model-generated density time series;
- model-generated constituent/temperature outputs;
- summary statistics and residual metrics;
- publication figures generated from those outputs;
- configuration/units metadata;
- hashes and citations for external inputs.

## What must not be published automatically

The result router fails closed on:

- credentials, tokens and private keys;
- private absolute machine paths;
- private application source;
- external coefficient/resource payloads whose policy is external-only;
- raw third-party datasets unless their redistribution permission has been
  explicitly established;
- Fortran/source archives from restricted authorities;
- temporary/cache/build artefacts.

## Required provenance

`RUN_MANIFEST.json` records the repository and Git identity used for the
calculation, the paper/execution identifiers, artefact SHA-256 values, and
source input provenance supplied by the publishing workflow.

This makes a result auditable without requiring redistribution of every raw
input used to produce it.
