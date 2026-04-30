# P1-P1 formulation

The stress tensor can be decomposed in a deviatoric and a sphecial part as follows

$$
\begin{equation}
    \pmb{\sigma} = \tilde{\pmb{\sigma}} + p \mathbf{I}
\end{equation}
$$

where $p = \frac{1}{3} \mathrm{tr}(\pmb{\sigma})$ is referred to as the mean stress, and $\tilde{\pmb{\sigma}}$ denotes the deviatoric stress tensor. 

Considering an initial stress $\pmb{\sigma}_0$, the above decomposition becomes

$$
\begin{equation}
    \pmb{\sigma} - \pmb{\sigma}_0 = \left( \tilde{\pmb{\sigma}} - \tilde{\pmb{\sigma}}_0 \right) + \left( p - p_0 \right) \mathbf{I}
    \label{eq:stress_0}
\end{equation}
$$


## Stress linearization
The deviatoric stress tensor relates to the deviatoric part of the elastic strain tensor through the shear modulus $G$, that is,

$$
\begin{equation}
    \tilde{\pmb{\sigma}}  = 2G \ \tilde{\pmb{\varepsilon}}_e = 2 G \left( \tilde{\pmb{\varepsilon}} - \tilde{\pmb{\varepsilon}}_{ne} \right),
    \label{eq:tilde_sigma_0}
\end{equation}
$$

where $\tilde{\pmb{\varepsilon}}$ and $\tilde{\pmb{\varepsilon}}_{ne}$ respectively denote the deviatoric parts of the total strain and non-elastic strain tensors.

The total volumetric strain and the deviatoric total strain tensor are given by

$$
\begin{align}
    &\varepsilon_v = \mathrm{tr} (\pmb{\varepsilon}) = \nabla \cdot \mathbf{u}
    \\
    &\tilde{\pmb{\varepsilon}} = \pmb{\varepsilon} - \frac{1}{3}\varepsilon_v\mathbf{I} = \frac{1}{2}\left( \nabla\mathbf{u} + \nabla\mathbf{u}^T \right) - \frac{1}{3} \left( \nabla \cdot \mathbf{u} \right) \mathbf{I},\label{eq:tilde_eps}
\end{align}
$$

The volumetric and the deviatoric parts of the non-elastic strain tensor can be respectively expressed as

$$
\begin{align}
    &\varepsilon_{ne,v} = \varepsilon_{ne,v}^k + \phi_2 \mathbf{F}^k : \delta \pmb{\sigma} - \phi_2 B_{ne,v}^k
    \label{eq:eps_nev}
    \\
    &\tilde{\pmb{\varepsilon}}_{ne} = \tilde{\pmb{\varepsilon}}_{ne}^k + \phi_2 \left( \tilde{\mathbb{G}}_{ne}^k : \delta \pmb{\sigma} - \tilde{\mathbf{B}}_{ne}^k \right)
    \label{eq:tilde_eps_ne_1}
\end{align}
$$

where $\delta \pmb{\sigma} = \pmb{\sigma} - \pmb{\sigma}^k$ (note that omitted superscript implies $k+1$), and

$$
\begin{align}
    &\varepsilon_{ne,v}^k = \mathrm{tr}(\pmb{\varepsilon}_{ne})
    \\
    &B_{ne,v}^k = \mathrm{tr}(\mathbf{B}_{ne}^k)
    \\
    &\mathbf{F}^k = 
    \begin{bmatrix}
        G_{ii11}^k & G_{ii12}^k & G_{ii13}^k \\
        G_{ii12}^k & G_{ii22}^k & G_{ii23}^k \\
        G_{ii13}^k & G_{ii23}^k & G_{ii33}^k
    \end{bmatrix}
    \label{eq:F_k}
    \\
    &\tilde{\pmb{\varepsilon}}_{ne}^k = \pmb{\varepsilon}_{ne}^k - \frac{1}{3} \varepsilon_{ne,v}^k\mathbf{I}
    \\
    &\tilde{\mathbf{B}}_{ne}^k = \mathbf{B}_{ne}^k - \frac{1}{3} B_{ne,v}^k\mathbf{I}
    \\
    &\tilde{\mathbb{G}}_{ne}^k = \mathbb{G}_{ne}^k - \frac{1}{3}\mathbf{I} \otimes \mathbf{F}^k
\end{align}
$$

> **Note:** See [1] for further details on the derivation of the above equations.

Substituting Eq. $\eqref{eq:tilde_eps_ne_1}$ into Eq. $\eqref{eq:tilde_sigma_0}$ leads to

$$
\begin{equation}
    \tilde{\pmb{\sigma}} = 2 G \left[ \tilde{\pmb{\varepsilon}} - \tilde{\pmb{\varepsilon}}_{ne}^k - \phi_2 \left( \tilde{\mathbb{G}}_{ne}^k : \delta \pmb{\sigma} - \tilde{\mathbf{B}}_{ne}^k \right) \right],
    \label{eq:tilde_sigma_1}
\end{equation}
$$

Finally, substituting Eq. $\eqref{eq:tilde_sigma_1}$ into Eq. $\eqref{eq:stress_0}$ and solving for $\pmb{\sigma}$ results in the following expression for the linearized stress tensor

$$
\begin{equation}
    \pmb{\sigma} = \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}}_0 + \tilde{\pmb{\varepsilon}} + \frac{1}{2G} p \mathbf{I} - \tilde{\pmb{\varepsilon}}_\text{rhs}^k \right)
\end{equation}
$$

where $\tilde{\pmb{\varepsilon}}_0 = \frac{1}{2G} \tilde{\sigma}_0$, and

$$
\begin{equation}
    \tilde{\mathbb{C}}_T^k = \left( \frac{1}{2G} \mathbb{I} + \phi_2 \tilde{\mathbb{G}}_{ne}^k \right)^{-1}
\end{equation}
$$



## Linearized momentum balance equation

$$
\begin{equation}
    \nabla \cdot \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}} + \frac{1}{2G} p \mathbf{I} \right)
    =
    \mathbf{b} + 
    \nabla \cdot \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}}_\text{rhs}^k - \tilde{\pmb{\varepsilon}}_0 \right)
\end{equation}
$$


## Linearized mass balance (mean stress) equation

$$
\begin{equation}
    \frac{1}{K}p = \pmb{\varepsilon}_{v} - \varepsilon_{ne,v}^k - \phi_2 \mathbf{F}^k : \delta \pmb{\sigma} + \phi_2 B_{ne,v}^k
    \label{eq:mean_stress_1}
\end{equation}
$$

We remind that $\delta \pmb{\sigma} = \delta \tilde{\pmb{\sigma}} + \delta p \mathbf{I}$. Now we consider that the deviatoric stress is freezed between two consecutive iterations, that is, $\delta \tilde{\pmb{\sigma}} = 0$, implying that $\delta \pmb{\sigma} = \delta p \mathbf{I}$. In this manner, Eq. $\eqref{eq:mean_stress_1}$ can be solved for $p$, resulting in

$$
\begin{equation}
    \frac{1}{K}p = \pmb{\varepsilon}_{v} - \varepsilon_{ne,v}^k - \phi_2 F_v^k \delta p + \phi_2 B_{ne,v}^k
    \label{eq:mean_stress_2}
\end{equation}
$$


## References
- [1] Honório, H.T., Franceschini, A., Ferronaro, M., Hajibeygi, H. *Salt cavern simulations with a stabilized mixed finite element formulation for low-order tetrahedral elements*. Computational Methods in Applied Mechanics and Engineering, (under review), 2026.