# Linear momentum

Let us consider a time domain $\mathcal{T} \in [t_0, t_f]$ and a spatial domain $\Omega \in \mathbb{R}^3$ bounded by a closed surface $\Gamma$ outward oriented by a normal vector $\hat{\mathbf{n}}$. The bounding surface can be further split into $\Gamma^u$ (for Dirichlet boundary conditions) and $\Gamma^\sigma$ (for Neumann boundary conditions), such that $\Gamma = \Gamma^u \cup \Gamma^\sigma$ and $\Gamma^u \cap \Gamma^\sigma = \emptyset$. 

The linear momentum balance equation for quasi-static loads can be written as

$$
\nabla \cdot \pmb{\sigma} = \rho \mathbf{g} \quad \forall \hspace{2mm} (\mathbf{x} \times t) \in (\Omega \times \mathcal{T}),
\label{eq:mom_0}
$$

where,

- $\pmb{\sigma}$: stress tensor $[\text{Pa}]$;
- $\rho$: density $[\text{kg}/\text{m}^3]$;
- $\mathbf{g}$: gravity acceleration vector $[\text{m}/\text{s}^2]$.

Equation $\eqref{eq:mom_0}$ is subjected to the following boundary and initial conditions:

$$
\begin{align}
    \mathbf{u}(\mathbf{x}, t) = \bar{\mathbf{u}}(\mathbf{x}, t) \quad &\forall \hspace{2mm} (\mathbf{x} \times t) \in (\Gamma^u \times \mathcal{T})  \nonumber
    \\
    \pmb{\sigma}(\mathbf{x}, t) \cdot \hat{\mathbf{n}} = \bar{\mathbf{t}}(\mathbf{x}, t) \quad &\forall \hspace{2mm} (\mathbf{x} \times t) \in (\Gamma^\sigma \times \mathcal{T})  \nonumber
    \\
    \pmb{\sigma}(\mathbf{x}, t_0) = \pmb{\sigma}_0(\mathbf{x}) \quad &\forall \hspace{2mm} \mathbf{x} \in \Omega  \nonumber
\end{align}
$$

where $\bar{\mathbf{u}}(\mathbf{x}, t)$ and $\bar{\mathbf{t}}(\mathbf{x}, t)$ are the displacement and traction vector functions prescribed at $\Gamma^u$ and $\Gamma^\sigma$, respectively, and $\pmb{\sigma}_0$ is the initial stress tensor field.

In Eq. $\eqref{eq:mom_0}$, the stress is calculated by Hooke's law, that is,

$$
\pmb{\sigma} = \pmb{\sigma}_0 + \mathbb{C}_e : \pmb{\varepsilon}_e
$$

where $\pmb{\varepsilon}_e$ is the elastic strain tensor, and $\mathbb{C}_e$ is the rank 4 elastic tensor (yellow spring in Fig. {fig_constitutive_model_0}). However, most constitutive models for geomaterials, especially salt rocks, comprise elastic, viscoelastic (i.e., time-dependent elastic), viscoplastic (i.e., time-dependent inelastic), and thermal strains.

Additionally, under small strain assumption, the kinematic relation between the total strain $\pmb{\varepsilon}$ and the displacement vector $\mathbf{u}$ is

$$
\pmb{\varepsilon} = \frac{1}{2} \left( \nabla \mathbf{u} + \nabla \mathbf{u}^T \right).
$$

Finally, the elastic strain $\pmb{\varepsilon}_e$ relates to the total strain $\pmb{\varepsilon}$ by the additive decomposition (valid under small strain assumption), that is

$$
\pmb{\varepsilon} = \pmb{\varepsilon}_{e} + \pmb{\varepsilon}_{ne} + \pmb{\varepsilon}_{th} \quad \rightarrow \quad \pmb{\varepsilon}_{e} = \pmb{\varepsilon} - \pmb{\varepsilon}_{ne} - \pmb{\varepsilon}_{th},
$$

where $\pmb{\varepsilon}_{ne}$ comprises all the time-independent elastic and inelastic strains. It can be represented as

<a id="eq-eps-i"></a>

$$
\pmb{\varepsilon}_{ne} = \sum_{i=1}^{N_{ne}} \pmb{\varepsilon}_{i},
$$

with $N_{ne}$ denoting the number of non-elastic elements included in the constitutive model. In this manner, the stress tensor can be expressed as

<a id="eq-stress-1"></a>

$$
\pmb{\sigma} = \pmb{\sigma}_0 + \mathbb{C}_e : \left( \pmb{\varepsilon} - \pmb{\varepsilon}_{ne} - \pmb{\varepsilon}_{th}\right).
\label{eq:stress_1}
$$

Expressions for the non-elastic strains $\pmb{\varepsilon}_{ne}$ depend on the constitutive model adopted. They are described in the Constitutive Models section. Additionally, these non-elastic strain often depend on the stress itself, which makes the problem to be non-linear. The stress linearization process is described in sections [Momentum P1]((../numerical_formulations/model_momentum_P1.md#sec-linearization-P1)) and [Momentum P1-P1](../numerical_formulations/model_momentum_P1.md#sec-linearization-P1P1).