# NRLMSIS 2.1 — Worked Workflow

## Example A: single atmospheric point

Assume a spacecraft state is known at a UTC epoch.

### 1. Establish the geometry

Record the epoch, altitude, latitude and longitude using the conventions required by the public interface. If local solar time is an explicit input, compute or provide it using a provenance-controlled method.

### 2. Establish space-weather inputs

Provide the model-specific solar and geomagnetic quantities. Preserve delayed, averaged and history semantics; do not replace them with convenient but scientifically different values.

### 3. Evaluate the public interface

Use the repository's public entry point rather than reaching into private/internal functions.

### 4. Validate outputs

At minimum check:

- finite temperature where returned;
- finite, positive total mass density in the physical model domain;
- constituent densities are finite/non-negative where returned;
- output units match the interface documentation.

### 5. Record provenance

Record the Git commit, input epoch/position, driver values/source, and any external coefficient/resource identity.

## Example B: orbital time series

Given authoritative trajectory samples

$$
(t_1,\mathbf{r}_1),\ldots,(t_N,\mathbf{r}_N),
$$

evaluate the model independently at every sample to obtain

$$
\rho_1,\ldots,\rho_N.
$$

Then aggregate only after pointwise evaluation. For a daily mean,

$$
\bar{\rho}_d =
\frac{1}{N_d}
\sum_{k=1}^{N_d}\rho_k.
$$

## Example C: interpretation for drag

If a downstream drag analysis uses

$$
a_D =
\frac{1}{2}
\rho
C_D
\frac{A}{m}
v_{\mathrm{rel}}^2,
$$

the atmospheric repository supplies $\rho$. The spacecraft ballistic parameters $C_D$, $A/m$ and the relative-flow velocity are external inputs and should not be attributed to the atmospheric model.
