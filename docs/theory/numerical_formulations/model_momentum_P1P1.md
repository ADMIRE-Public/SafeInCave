# P1-P1 formulation

<!-- The stress tensor can be decomposed in a deviatoric and a sphecial part as follows -->

<!-- $$
\begin{equation}
    \pmb{\sigma} = \tilde{\pmb{\sigma}} + p \mathbf{I}
    \label{eq:stress_0}
\end{equation}
$$ -->

The stress tensor can be represented in terms of the elastic strain tensor as

$$
\begin{equation}
    \pmb{\sigma} = \pmb{\sigma}_0 + 2G \pmb{\varepsilon}_e + \lambda \varepsilon_{e,v} \mathbf{I}.
\end{equation}
$$

where $G$ and $\lambda$ denote the shear modulus and first Lamè's parameter, respectively. Decomposing the elastic strain tensor in its deviatoric $\tilde{\pmb{\varepsilon}}_e$ and spherical parts, gives

$$
\begin{equation}
    \pmb{\sigma} = \pmb{\sigma}_0 + 2G \left( \tilde{\pmb{\varepsilon}}_e + \frac{1}{3} \varepsilon_{e,v} \right) + \lambda \varepsilon_{e,v} \mathbf{I}
\end{equation}
$$

Recognizing that the elastic bulk modulus is given by $K = \frac{2G}{3} + \lambda$, and the mean stress relates to the volumetric elastic strain by $p = K \varepsilon_{e,v}$, it follows that,

$$
\begin{equation}
    \pmb{\sigma} = \pmb{\sigma}_0 + 2G \tilde{\pmb{\varepsilon}}_e + p \mathbf{I}
\end{equation}
$$

Moreover, we have defined that $\tilde{\pmb{\varepsilon}}_e = \tilde{\pmb{\varepsilon}} - \tilde{\pmb{\varepsilon}}_{ne}$, therefore, the stress tensor can be represented as

$$
\begin{equation}
    \pmb{\sigma} = \pmb{\sigma}_0 + 2G \left( \tilde{\pmb{\varepsilon}} - \tilde{\pmb{\varepsilon}}_{ne} \right) + p \mathbf{I}
    \label{eq:stress_0}
\end{equation}
$$

where $\tilde{\pmb{\varepsilon}}_{ne}$ is the deviatoric part of the non-elastic strain tensor. The term $\tilde{\pmb{\varepsilon}}_{ne}$ is a function of the stress itself (among other parameters). As result, Eq. $\eqref{eq:stress_0}$ becomes non-linear. The linearization of Eq. $\eqref{eq:stress_0}$ is described below.

<!-- where $p = \frac{1}{3} \mathrm{tr}(\pmb{\sigma})$ is referred to as the mean stress, and $\tilde{\pmb{\sigma}}$ denotes the deviatoric stress tensor.  -->

<!-- Considering an initial stress $\pmb{\sigma}_0$, the above decomposition becomes

$$
\begin{equation}
    \pmb{\sigma} - \pmb{\sigma}_0 = \left( \tilde{\pmb{\sigma}} - \tilde{\pmb{\sigma}}_0 \right) + \left( p - p_0 \right) \mathbf{I}
    \label{eq:stress_0}
\end{equation}
$$ -->


<a id="sec-linearization-P1P1"></a>
## Stress linearization

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

??? note "Note"

    See [1] for further details on the derivation of the above equations.

Substituting Eq. $\eqref{eq:tilde_eps_ne_1}$ into Eq. $\eqref{eq:stress_0}$ and dividing by $2G$ leads to

$$
\begin{equation}
    \frac{1}{2G}\pmb{\sigma} = \frac{1}{2G}\pmb{\sigma}_0 + \tilde{\pmb{\varepsilon}} - \tilde{\pmb{\varepsilon}}_{ne}^k - \phi_2 \left( \tilde{\mathbb{G}}_{ne}^k : \delta \pmb{\sigma} - \tilde{\mathbf{B}}_{ne}^k \right) + p \mathbf{I},
\end{equation}
$$

Recalling that $\delta \pmb{\sigma} = \pmb{\sigma} - \pmb{\sigma}^k$, then

$$
\begin{equation}
    \left( \frac{1}{2G} \mathbb{I} + \phi_2 \tilde{\mathbb{G}}_{ne}^k \right) : \pmb{\sigma} = \tilde{\pmb{\varepsilon}} + \frac{1}{2G} p \mathbf{I} + \frac{1}{2G} \pmb{\sigma}_0 - \tilde{\pmb{\varepsilon}}_\text{rhs}^k
\end{equation}
$$

Finally, the linearized stress tensor can be expressed as

$$
\begin{equation}
    \pmb{\sigma} = \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}} + \frac{1}{2G} p \mathbf{I} + \frac{1}{2G} \pmb{\sigma}_0 - \tilde{\pmb{\varepsilon}}_\text{rhs}^k \right).
    \label{eq:stress_linearized}
\end{equation}
$$



## Linearized momentum balance equation

The momentum balance equation is written as

$$
\nabla \cdot \pmb{\sigma} = \mathbf{b}.
\label{eq:mom_0}
$$

Substituting Eq. $\eqref{eq:stress_linearized}$ into Eq. $\eqref{eq:mom_0}$, leads to the linearized momemtum balance equation, which can be expressed as

$$
\begin{equation}
    \nabla \cdot \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}} + \frac{1}{2G} p \mathbf{I} \right)
    =
    \mathbf{b} + 
    \nabla \cdot \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}}_\text{rhs}^k - \frac{1}{2G} \pmb{\sigma}_0 \right),
\end{equation}
$$

where the total volumetric strain and the deviatoric total strain tensors are given by

$$
\begin{align}
    &\varepsilon_v = \mathrm{tr} (\pmb{\varepsilon}) = \nabla \cdot \mathbf{u},
    \\
    &\tilde{\pmb{\varepsilon}} = \pmb{\varepsilon} - \frac{1}{3}\varepsilon_v\mathbf{I} = \frac{1}{2}\left( \nabla\mathbf{u} + \nabla\mathbf{u}^T \right) - \frac{1}{3} \left( \nabla \cdot \mathbf{u} \right) \mathbf{I}.
    \label{eq:tilde_eps}
\end{align}
$$


## Linearized mass balance (mean stress) equation

In solid mechanics, the mass balance equation relates the mean stress $p$ to the elastic volumetric strain $\varepsilon_{e,v}$, that is,

$$
p = K \varepsilon_{e,v} = K \left( \varepsilon_{v} - \varepsilon_{th,v} - \varepsilon_{ne,v} \right)
\label{eq:eq:mean_stress_0}
$$

where $K$ is the elastic bulk modulus, $\pmb{\varepsilon}_{v}$ is the total volumetric strain, and $\varepsilon_{ne,v}$ is the non-elastic volumetric strain.


Substituting Eq. $\eqref{eq:eps_nev}$ into Eq. $\eqref{eq:eq:mean_stress_0}$,

$$
p = K \left( \varepsilon_{v} - \varepsilon_{th,v} - \varepsilon_{ne,v}^k - \phi_2 \mathbf{F}^k : \delta \pmb{\sigma} + \phi_2 B_{ne,v}^k \right)
\label{eq:mean_stress_1}
$$

We remind that $\delta \pmb{\sigma} = \delta \tilde{\pmb{\sigma}} + \delta p \mathbf{I}$. We now consider that the deviatoric stress is freezed between two consecutive iterations, that is, $\delta \tilde{\pmb{\sigma}} = 0$, implying that $\delta \pmb{\sigma} = \delta p \mathbf{I}$. In this manner, Eq. $\eqref{eq:mean_stress_1}$ can be solved for $p$, resulting in

$$
\begin{equation}
    K^{-1}p = \varepsilon_{v} - \varepsilon_{th,v} - \varepsilon_{ne,v}^k - \phi_2 F_v^k \delta p + \phi_2 B_{ne,v}^k
    \label{eq:mean_stress_2}
\end{equation}
$$

where $F_v^k = \mathbf{F}^k : \mathbf{I} = \mathrm{tr}(\mathbf{F}^k)$. Again, recalling that $\delta p = p - p^k$, after some manipulation, the linearized mean stress equation can be expressed as

$$
\begin{equation}
    K_T^k p - \varepsilon_{v} = b_p^k - \varepsilon_{ne,v}^k - \varepsilon_{th,v}
    \label{eq:mean_stress_3}
\end{equation}
$$

where $b_p^k = \phi_2 \left(F_v^k p^k + B_{ne,v}^k\right)$ and $K_T^k = K^{-1} + \phi_2 F_v^k$.


<!-- ## Weak formulation
For the mixed formulation, the displacement field $\mathbf{u}$ is approximated in a first-order Sobolev space $H^1(\Omega)$, while mean stress field $p$ is sufficient to be defined in $L^2(\Omega)$. In this manner, the continuous trial function spaces are defined as

$$
\begin{align}
    &\mathcal{U} = \lbrace \mathbf{u} : \Omega \rightarrow \mathbb{R}^3 \hspace{1mm} | \hspace{1mm} \mathbf{u} \in [H^1(\Omega)]^3, \mathbf{u} = \bar{\mathbf{u}} \hspace{1mm} \text{on} \hspace{1mm} \Gamma^u \rbrace, \nonumber
    \\
    &\mathcal{P} = \lbrace p : \Omega \rightarrow \mathbb{R}\hspace{1mm} | \hspace{1mm} p \in L^2(\Omega) \rbrace, \nonumber
\end{align}
$$

Since no essential boundary conditions are imposed on the mean stress $p$, an unconstrained $L^2(\Omega)$ space is sufficient for both the trial and test functions. The continuous test (weighting) function space for the displacement field is

$$
\mathcal{U}^0 = \lbrace \mathbf{w} : \Omega \rightarrow \mathbb{R}^3 \hspace{1mm} | \hspace{1mm} \mathbf{w} \in [H^1(\Omega)]^3, \mathbf{w} = \mathbf{0} \hspace{1mm} \text{on} \hspace{1mm} \Gamma^u \rbrace, \nonumber
$$

$$
\begin{split}
	\int_\Omega \pmb{\varepsilon}(\mathbf{w}) : \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}} + \frac{1}{2G}p\mathbf{I} \right) \mathrm{d}\Omega
	=
	\int_\Omega \rho\mathbf{g} \cdot \mathbf{w} \mathrm{d}\Omega
	+
	\int_\Gamma \mathbf{t} \cdot \mathbf{w} \mathrm{d} \Gamma
	+
	\\
	+ \int_\Omega \pmb{\varepsilon}(\mathbf{w}) : \tilde{\mathbb{C}}_T^k : \left( \tilde{\pmb{\varepsilon}}_\text{rhs}^k - \tilde{\pmb{\varepsilon}}_0 \right) \mathrm{d}\Omega,
	\quad \forall \hspace{1mm} \mathbf{w} \in \mathcal{U}^0.
\end{split}
$$

where

$$
\tilde{\pmb{\varepsilon}} = \pmb{\varepsilon} - \frac{1}{3} \varepsilon_v \mathbf{I} = \frac{1}{2} \left( \nabla \mathbf{u} + \nabla \mathbf{u}^T \right) - \frac{1}{3} \left(\nabla \cdot \mathbf{u}\right) \mathbf{I}
$$

and the unknowns are $\mathbf{u}$ and $p$. -->

## References
- [1] Honório, H.T., Franceschini, A., Ferronaro, M., Hajibeygi, H. *Salt cavern simulations with a stabilized mixed finite element formulation for low-order tetrahedral elements*. Computational Methods in Applied Mechanics and Engineering, 2026.