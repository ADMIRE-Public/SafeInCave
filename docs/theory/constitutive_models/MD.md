# Munson-Dawson Constitutive Model

>**Note**: Mention here that the text below was largely based on Gijs' Master thesis.

The Munson-Dawson model describes salt creep through a unified inelastic strain-rate law governed by internal state variables.

## Steady-state creep
The steady-state creep rate follows a thermally activated power-law relation,

$$
\dot{\varepsilon}_{ss}
=
A
\exp\!\left(-\frac{Q}{RT}\right)
q^n,
$$

where $q$ is the Von Mises stress.

## Unified inelastic strain rate
The inelastic strain-rate tensor is defined as

$$
\dot{\varepsilon}_{ij}
=
F(\zeta,q)
\,
\dot{\varepsilon}_{ss}
\,
\frac{\partial q}{\partial \sigma_{ij}},
$$

with

$$
\frac{\partial q}{\partial \sigma_{ij}}
=
\frac{3}{2}
\frac{s_{ij}}{q}.
$$

This ensures that the strain-rate tensor remains coaxial with the deviatoric stress tensor.

## Transient multiplier and internal variable
The transient multiplier is defined piecewise as

$$
F =
\begin{cases}
    \exp\!\left[
    \Delta \left(1-\frac{\zeta}{\varepsilon_t^*}\right)^2
    \right],
    & \zeta < \varepsilon_t^*,
    \\
    1,
    & \zeta = \varepsilon_t^*,
    \\
    \exp\!\left[
    -\delta \left(1-\frac{\zeta}{\varepsilon_t^*}\right)^2
    \right],
    & \zeta > \varepsilon_t^*,
\end{cases}
$$

where $\delta$ is a recovery parameter and

$$
\Delta = \alpha_w + \beta_w \log_{10}\!\left(\frac{q}{\mu}\right)
$$

is a stress-dependent hardening parameter, with $\alpha_w$ and $\beta_w$ material constants.

The internal variable $\zeta$ represents the accumulated transient strain and evolves
according to

$$
\dot{\zeta}
=
(F-1)\dot{\varepsilon}_{ss}.
$$

This evolution is self-regulating. When the stress state changes, $\zeta$ is no longer at
its equilibrium value and $F \neq 1$, producing a transient creep response. Since $F > 1$
implies $\dot{\zeta} > 0$, the internal variable grows toward its equilibrium
$\varepsilon_t^*$, which progressively reduces $F$ toward unity. At equilibrium
($\zeta = \varepsilon_t^*$), $F = 1$ and $\dot{\zeta} = 0$, so the model returns to
steady-state creep. Conversely, when $\zeta > \varepsilon_t^*$ (e.g.\ after a stress
reduction), $F < 1$ and $\dot{\zeta} < 0$, so that $\zeta$ decreases back toward the new
equilibrium.

The saturation strain is

$$
\varepsilon_t^*
=
K_0 e^{cT}
\left(
\frac{q}{\mu}
\right)^m,
$$

where $K_0$, $c$, and $m$ are material parameters and $\mu$ is the shear modulus. Because
$\varepsilon_t^*$ depends on the current stress level, a change in loading conditions shifts
the equilibrium and triggers a new transient phase.