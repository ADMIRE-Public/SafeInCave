# Heat diffusion
Let us define the following trial function space for approximating temperature as

$$
\mathcal{S} = \lbrace T : \Omega \rightarrow \mathbb{R}\hspace{1mm} | \hspace{1mm} T \in H^1, T = \bar{T} \hspace{1mm} \text{on} \hspace{1mm} \Gamma_d^T \rbrace,
\label{eq:trial_T}
$$


and the test function space as

$$
\mathcal{S}_0 = \lbrace v : \Omega \rightarrow \mathbb{R}\hspace{1mm} | \hspace{1mm} v \in H^1, v = 0 \hspace{1mm} \text{on} \hspace{1mm} \Gamma_d^T \rbrace.
\label{eq:test_T}
$$

The weak form of the heat diffusion equation reads

$$
\begin{align}
    \int_\Omega \left( \rho c \frac{\partial T}{\partial t} + k\nabla T \cdot \nabla v \right)\mathrm{d}\Omega 
    + \int_{\Gamma^q} q'' v \mathrm{d}\Gamma
    + \int_{\Gamma^h} h\left( T - T_\infty \right) v \mathrm{d}\Gamma
    = 0.
    \label{eq:weak_heat_0}
\end{align}
$$


Integrating in time between $t$ and $t+\Delta t$, using the fully-implicit (i.e., backward Euler) scheme to evaluate the time  integrals, and rearranging the terms, Eq. $\eqref{eq:weak_heat_0}$ becomes

$$
\begin{align}
    \int_\Omega \left(  \frac{\rho c}{\Delta t}T + k\nabla T \cdot \nabla v \right)\mathrm{d}\Omega 
    + \int_{\Gamma^h} h T v \mathrm{d}\Gamma
    = 
    \int_\Omega  \frac{\rho c}{\Delta t} T^t \mathrm{d}\Omega 
    + \int_{\Gamma^h} h T_\infty v \mathrm{d}\Gamma
    - \int_{\Gamma^q} q'' v \mathrm{d}\Gamma
    . \nonumber
    \label{eq:weak_heat_1}
\end{align}
$$

where $T^t$ refers to the temperature evaluated at the previous time step $t$, while the temperature evaluated at the current time step $t+\Delta t$ carries no superscript to avoid heavy notation ($T$). 