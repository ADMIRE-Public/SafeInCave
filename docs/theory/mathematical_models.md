# Governing Equations

This section presents the governing equations considered in the SafeInCave simulator, namely, the heat diffusion equation, the linear momentum balance equations, and the thermodynamic model for brine (used for cavern abandonment simulations).

## Heat diffusion

Let us consider a time domain $\mathcal{T} \in [t_0, t_f]$ and a spatial domain $\Omega \in \mathbb{R}^3$ bounded by a closed surface $\Gamma$ outward oriented by a normal vector $\hat{\mathbf{n}}$. Consider the surface $\Gamma$ can be further split into $\Gamma^T$, $\Gamma^q$, and $\Gamma^h$, such that $\Gamma = \Gamma^T \cup \Gamma^q \cup \Gamma^h$, and $\Gamma^T \cap \Gamma^q = \Gamma^T \cap \Gamma^h = \Gamma^h \cap \Gamma^q = \emptyset$. The heat diffusion equation without heat generation can be expressed as

$$
\begin{equation}
    \rho c \frac{\partial T}{\partial t} - \nabla \cdot \left( k\nabla T \right) = 0 \quad \forall \hspace{2mm} 
    (\mathbf{x} \times t) \in (\Omega \times \mathcal{T})
    \label{eq:heat_0}
\end{equation}
$$

where,

- $T$: temperature $[\text{K}]$;
- $c$: specific heat capacity $[\text{J}/\text{kg}/\text{K}]$;
- $\rho$: density $[\text{kg}/\text{m}^3]$;
- $k$: thermal conductivity $[\text{W}/\text{m}^3]$.

> **_NOTE:_** The thermal conductivity $k$ is strictly considered to be a scalar here, although it can vary in space (as any other material property).

Equation $\eqref{eq:heat_0}$ is subjected to the following boundary and initial conditions:

$$
\begin{align}
    T(\mathbf{x}, t) = \bar{T}(\mathbf{x}, t) \quad &\forall \hspace{2mm} (\mathbf{x} \times t) \in (\Gamma^T \times \mathcal{T}) \nonumber
    \\
    -k\nabla T(\mathbf{x}, t) \cdot \hat{\mathbf{n}} = \bar{q}''(\mathbf{x}, t) \quad &\forall \hspace{2mm} (\mathbf{x} \times t) 
    \in (\Gamma^q \times \mathcal{T})  \nonumber
    \\
    -k\nabla T(\mathbf{x}, t) \cdot \hat{\mathbf{n}} = h\left( T - T_\infty \right) \quad &\forall \hspace{2mm} (\mathbf{x} \times t) 
    \in (\Gamma^h \times \mathcal{T})  \nonumber
    \\
    T(\mathbf{x}, t_0) = T_0(\mathbf{x}) \quad &\forall \hspace{2mm} \mathbf{x} \in \Omega  \nonumber
\end{align}
$$

where $\bar{T}(\mathbf{x}, t)$ and $\bar{q}''(\mathbf{x}, t)$ are the temperature and heat flux functions prescribed at $\Gamma^T$ and $\Gamma^q$, respectively. Additionally, $h$ is the convective heat transfer coefficient [$\text{W}/\text{m}^2/\text{K}$], and $T_\infty$ is the far field temperature (usually the gas/brine temperature).


## Linear momentum

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

where $\pmb{\varepsilon}_e$ is the elastic strain tensor, and $\mathbb{C}_e$ is the rank 4 elastic tensor (yellow spring in Fig. {fig_constitutive_model_0}). However, most constitutive models for geomaterials, especially salt rocks, comprise elastic, viscoelastic (i.e., time-dependent elastic), and viscoplastic (i.e., time-dependent inelastic) deformations.

>**_NOTE:_** The term non-elastic deformation includes all types of deformation that are not instantaneously elastic, that is, viscoelastic (time dependent elastic) and inelastic (viscoplastic, plastic, creep, etc) deformations.

Small strain assumption is adopted, so that the additive decomposition holds for the total strain tensor, that is

$$
\pmb{\varepsilon} = \pmb{\varepsilon}_{e} + \pmb{\varepsilon}_{ne} + \pmb{\varepsilon}_{th} \quad \rightarrow \quad \pmb{\varepsilon}_{e} = \pmb{\varepsilon} - \pmb{\varepsilon}_{ne} - \pmb{\varepsilon}_{th},
$$

where the kinematic relation for small strains provides

$$
\pmb{\varepsilon} = \frac{1}{2} \left( \nabla \mathbf{u} + \nabla \mathbf{u}^T \right),
$$

and the non-elastic strains are given by

$$
\pmb{\varepsilon}_{ne} = \sum_{i=1}^{N_{ne}} \pmb{\varepsilon}_{i},
$$

with $N_{ne}$ denoting the number of non-elastic elements included in the constitutive model. In this manner, the stress tensor can be expressed as

$$
\pmb{\sigma} = \pmb{\sigma}_0 + \mathbb{C}_e : \left( \pmb{\varepsilon} - \pmb{\varepsilon}_{ne} - \pmb{\varepsilon}_{th}\right).
\label{eq:stress_1}
$$

From Section [Constitutive models](constitutive_models.md#constitutive-model), it was shown that the non-elastic elements have expressions for the strain rates (i.e., $\dot{\pmb{\varepsilon}}_{ne}$), not for the strains $\pmb{\varepsilon}_{ne}$, as required by Eq. $\eqref{eq:stress_1}$. This will be addressed later in this section.
