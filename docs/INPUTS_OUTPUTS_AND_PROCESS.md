# NRLMSIS 2.1 — Inputs, Outputs and Computation Process

## Required input categories

- UTC epoch / day of year and time of day
- geodetic altitude, latitude and longitude
- local apparent solar time
- previous-day F10.7
- 81-day mean F10.7
- geomagnetic activity input/history required by the accepted interface

The exact argument names, units, valid ranges and optional modes are defined by the repository's public Python interface and tests.

## Outputs

- neutral temperature and exospheric temperature
- major neutral-species number densities
- nitric-oxide number density in the NRLMSIS 2.1 domain
- total neutral mass density

## Computation sequence

1. Validate the requested epoch, position and model domain.
2. Resolve or receive the required solar and geomagnetic drivers.
3. Convert the state to the conventions expected by the scientific core.
4. Evaluate model-specific temporal, spatial and environmental basis functions.
5. Evaluate temperature and constituent/reference-density structure.
6. Construct total neutral mass density.
7. Convert outputs to the public interface units.
8. Preserve provenance for the driver source, model commit and any external resource.

## Time-series use

For a trajectory with states $\{(t_k,\mathbf{r}_k)\}_{k=1}^N$, evaluate the model at each authoritative epoch/position:

```math
\rho_k =
\mathcal{M}_\rho
\left(t_k,\mathbf{r}_k,\mathbf{s}_k\right).
```

A daily product is then derived from the valid pointwise evaluations; it is not a different atmospheric model.

## Arithmetic daily mean

For a UTC day with $N_d$ valid samples,

```math
\bar{\rho}_d =
\frac{1}{N_d}
\sum_{k=1}^{N_d} \rho_k.
```

Coverage metadata should accompany the mean whenever the underlying trajectory is incomplete.
