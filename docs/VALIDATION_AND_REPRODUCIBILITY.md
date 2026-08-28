# NRLMSIS 2.1 — Validation and Reproducibility

## Validation principle

Three identities should be kept separate:

1. **scientific identity** — whether the public numerical core is the accepted model implementation;
2. **interface identity** — whether public inputs/outputs preserve the accepted semantics and units;
3. **application validation** — comparison against an external density/temperature reference.

A good external comparison does not replace implementation-equivalence evidence, and an implementation-equivalence result does not imply perfect agreement with nature.

## Repository reproducibility record

A reproducible calculation should capture:

- repository URL and Git commit;
- Python/runtime versions;
- exact input epochs and positions;
- solar/geomagnetic driver source and values;
- external coefficient/resource identity where applicable;
- public interface name and units;
- result CSV/figure hashes where publication products are produced.

## Public tests

The repository's existing tests remain part of the acceptance boundary. This documentation enhancement must not change the scientific source files or make an existing public test fail.

## External comparisons

When comparing with an external reference, report the reference provenance and sampling coverage. A useful signed relative difference is

```math
\delta_\rho =
\frac{\rho_{\mathrm{model}}-\rho_{\mathrm{ref}}}
{\rho_{\mathrm{ref}}}.
```

The corresponding percentage is

```math
100\,\delta_\rho.
```

These are comparison metrics, not automatic declarations that either source is defective.

## Coverage-aware daily validation

A complete-day statistic and an available-sample statistic must not be silently conflated. Store $N_d$ and, where a nominal cadence is known, the coverage fraction

```math
C_d = \frac{N_d}{N_{\mathrm{expected},d}}.
```

## Accepted Elara X lineage

This repository was published only after controlled scientific/provenance acceptance in the Elara X atmospheric-model programme. The enhancement phase changes documentation and GitHub presentation only; accepted scientific code remains immutable.
