# Heat diffusion

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

