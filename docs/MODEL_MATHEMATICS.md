# NRLMSIS 2.1 — Governing Mathematics

## Purpose

NRLMSIS 2.0 thermosphere/whole-atmosphere formulation with the NRLMSIS 2.1 nitric-oxide extension.

The equations below explain the physical and mathematical structure of the model. They are intentionally separated from the accepted implementation source: documentation must not silently become a second, divergent implementation.

## Empirical state dependence

A thermospheric empirical model can be represented schematically as a mapping

$$
\mathcal{M}:
(t,\mathbf{r},\mathbf{s})
\longrightarrow
(T,\rho,\mathbf{n}),
$$

where $t$ is epoch, $\mathbf{r}$ is position, $\mathbf{s}$ contains the required space-weather drivers, $T$ is temperature, $\rho$ is total mass density and $\mathbf{n}$ denotes constituent densities where provided.


## Thermospheric temperature structure

A useful representation of the upper-atmosphere temperature branch is the Bates-type form

$$
T(z) = T_\infty - \left(T_\infty - T_\ell\right)
\exp\!\left[-s\left(z-z_\ell\right)\right],
$$

where $T_\infty$ is exospheric temperature, $T_\ell$ is the temperature at a lower reference level, and $s$ controls the vertical temperature gradient.

The precise accepted implementation contains the fitted model basis and transition logic. This equation is therefore explanatory rather than a replacement implementation.

## Hydrostatic/diffusive structure

For a constituent $i$, hydrostatic balance can be written as

$$
\frac{dp_i}{dz} = -\rho_i g,
$$

with

$$
p_i = n_i k_B T,
\qquad
\rho_i = m_i n_i.
$$

Above the mixed lower thermosphere, the constituent profiles approach species-dependent diffusive behaviour; lower down, the model transitions toward mixed-atmosphere behaviour using the formulation encoded in the accepted scientific implementation.

## Total density

Total mass density is obtained from the constituent number densities:

$$
\rho = \sum_i m_i n_i.
$$

The species set and any model-version-specific terms are defined by the accepted implementation.

## NRLMSIS 2.1 nitric oxide extension

NRLMSIS 2.1 retains the NRLMSIS 2.0 temperature and non-NO species formulation and adds a fitted nitric-oxide vertical profile. A schematic representation is

$$
n_{\mathrm{NO}}(z,\mathbf{x}) =
n_{\mathrm{NO},0}(z)
\exp\!\left(G_{\mathrm{NO}}(\mathbf{x})\right),
$$

where $\mathbf{x}$ represents the environmental variables used by the NO fit. The public documentation should be read together with the accepted source and the NRLMSIS 2.1 scientific paper; this schematic equation is not a substitute for the fitted coefficient implementation.


## Unit discipline

Density is normally exposed by the Elara X public interfaces in SI units of kg m$^{-3}$ unless the interface explicitly documents a model-native unit. Angles, altitude and space-weather indices must follow the repository interface contract. Do not infer units from a variable name alone.

## Scientific reference

Emmert et al. (2022), NRLMSIS 2.1: An Empirical Model of Nitric Oxide Incorporated Into MSIS, DOI 10.1029/2022JA030896.
