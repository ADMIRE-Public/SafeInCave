
# P1 formulation

## Stress linearization
Consider a non-elastic element $i$, with strain rate $\dot{\pmb{\varepsilon}}_i$. Using the $\theta-$method to integrate the strain rate from $t$ to $t+\Delta t$, gives

$$
\pmb{\varepsilon}_i = \pmb{\varepsilon}_i^t + \underbrace{\Delta t \theta}_{\phi_1} \dot{\pmb{\varepsilon}}_i^t + \underbrace{\Delta t (1 -\theta)}_{\phi_2} \dot{\pmb{\varepsilon}}_i.
\label{eq:eps_i_0}
$$

The strain rate $\dot{\pmb{\varepsilon}}_i$, evaluated at the current time level $t+\Delta t$, depends on the stress tensor $\pmb{\sigma}$, which is still an unknown. Therefore, the problem must be linearized and solved iteratively in each time step. For this reason, for every variable at $t+\Delta t$ we indicate the **previous** iteration level with the superscript $k$ (previous iteration), and to the **current** iteration with $k+1$. Therefore, Eq. $\eqref{eq:eps_i_0}$ is rewritten as

$$
\pmb{\varepsilon}_i^{k+1} = \pmb{\varepsilon}_i^t + \phi_1 \dot{\pmb{\varepsilon}}_i^t + \phi_2 \dot{\pmb{\varepsilon}}_i^{k+1}
\label{eq:eps_i_1}
$$

From Equations [(6)](mathematical_models.md#eq-eps-i) and [(7)](mathematical_models.md#eq-stress-1), it follows that the stress tensor at the current iteration $k+1$ can be expressed as

$$
\pmb{\sigma}^{k+1} = \pmb{\sigma}_0 + \mathbb{C}_e : \left( \pmb{\varepsilon}^{k+1} - \pmb{\varepsilon}_{th} - \pmb{\varepsilon}_{ne}^t - \phi_1 \dot{\pmb{\varepsilon}}_{ne}^t - \phi_2 \dot{\pmb{\varepsilon}}_{ne}^{k+1}\right).
\label{eq:stress_2}
$$

>**_NOTE:_** No superscript is used for the thermal strain $\pmb{\varepsilon}_{th}$ as it only depends on temperature, so the nonlinear iterations do not apply. But keep in mind it refers to the current time level $t+\Delta t$.

$$
\dot{\pmb{\varepsilon}}_i^{k+1} = \dot{\pmb{\varepsilon}}_i^{k} 
	+ \frac{\partial \dot{\pmb{\varepsilon}}_i}{\partial \pmb{\sigma}} : \delta \pmb{\sigma}
	+ \frac{\partial \dot{\pmb{\varepsilon}}_i}{\partial \omega_i} \delta \omega_i
\label{eq:eps_rate_i_0}
$$

where $\delta\pmb{\sigma} = \pmb{\sigma}^{k+1} - \pmb{\sigma}^k$ and $\delta\omega_i = \omega_i^{k+1} - \omega_i^k$.

>**_NOTE:_** The term $\frac{\partial \dot{\pmb{\varepsilon}}_i}{\partial \pmb{\sigma}}$ is a rank-4 tensor, while $\delta\pmb{\sigma}$ is arank-2, hence the double dot product between them, which results a rank-2 tensor. For further support on tensorial operations, check [here](https://youtu.be/w5KX3F_rdzU?si=QQLVBq1NcrvOiS32), and [here](https://youtu.be/JiN6jwp0RPk?si=K1Qhe3lAxJD4LI5w) for practical examples.

The increment of state variable $\delta\omega_i$ can be obtained by defining a residual equation based on the evolution equation of $\omega_i$ and using Newton-Raphson to drive the residual to zero. Considering the residual equation is of the form $r_i = r_i(\pmb{\sigma}, \omega_i)$, it follows that

$$
r_i^{k+1} = r_i^k + \frac{\partial r_i}{\partial \pmb{\sigma}} : \delta \pmb{\sigma} + \underbrace{\frac{\partial r_i}{\partial \omega_i}}_{h_i} \delta \omega_i = 0
	\quad \rightarrow \quad
	\delta \omega_i = - \frac{1}{h_i} \left( r_i^k + \frac{\partial r_i}{\partial \pmb{\sigma}} : \delta \pmb{\sigma} \right).
\label{eq:res_0}
$$

>**_NOTE:_** Currently, only the viscoplastic element uses a state variable (the hardening parameter, $\alpha$). For this case, the residual equation would read $r(\alpha) = \alpha - a_1 \left[ \left( a_1 / \alpha_0 \right)^{1/\eta} + \xi \right]^{-1}$.

Substituting Eq. $\eqref{eq:res_0}$ into Eq. $\eqref{eq:eps_rate_i_0}$ to eliminate $\delta\omega_i$ yields

$$
\dot{\pmb{\varepsilon}}_i^{k+1} = \dot{\pmb{\varepsilon}}_i^{k} 
	+ \underbrace{\left( \frac{\partial \dot{\pmb{\varepsilon}}_i}{\partial \pmb{\sigma}} - \frac{1}{h_i} \frac{\partial \dot{\pmb{\varepsilon}}_i}{\partial \omega_i} \frac{\partial r_i}{\partial \pmb{\sigma}} \right)}_{\mathbb{G}_i} : \delta \pmb{\sigma}
	- \underbrace{\frac{r_i^k}{h_i} \frac{\partial \dot{\pmb{\varepsilon}}_i}{\partial \omega_i}}_{\mathbf{B}_i}
	\quad \rightarrow \quad
	\dot{\pmb{\varepsilon}}_i^{k+1} = \dot{\pmb{\varepsilon}}_i^{k} + \mathbb{G}_i : \delta \pmb{\sigma} - \mathbf{B}_i.
\label{eq:eps_rate_i_1}
$$

Considering all non-elastic elements,

$$
\dot{\pmb{\varepsilon}}_{ne}^{k+1}
	=
	\dot{\pmb{\varepsilon}}_{ne}^{k} + \mathbb{G}_{ne} : \delta \pmb{\sigma} - \mathbf{B}_{ne},
\label{eq:eps_rate_ne_0}
$$

where $\mathbb{G}_{ne} = \sum_{i=1}^{N_{ne}} \mathbb{G}_i$ and $\mathbf{B}_{ne} = \sum_{i=1}^{N_{ne}} \mathbf{B}_i$.

Finally, substituting Eq. $\eqref{eq:eps_rate_ne_0}$ into Eq. $\eqref{eq:stress_2}$ leads to

$$
\pmb{\sigma}^{k+1} = \mathbb{C}_T : \left[ 
		\pmb{\varepsilon}_0 
		+ \pmb{\varepsilon}^{k+1} 
		- \pmb{\varepsilon}_{ne}^k
		+ \phi_2 \left( \mathbb{G}_{ne} : \pmb{\sigma}^k + \mathbf{B}_{ne} \right)
	\right]
\label{eq:stress_3}
$$

where the initial elastic strain (consistent with the imposed initial stress field), the non-elastic strain at previous iteration $k$, and the consistent tangent matrix are respectively given by

$$
\begin{align}
    &\pmb{\varepsilon}_0 = \mathbb{C}_e^{-1} : \pmb{\sigma}_0
    \\
    &\pmb{\varepsilon}_{ne}^k = \pmb{\varepsilon}_{ne}^t + \phi_1 \dot{\pmb{\varepsilon}}_{ne}^t + \phi_1 \dot{\pmb{\varepsilon}}_{ne}^k
    \\
    &\mathbb{C}_T = \left( \mathbb{C}_0^{-1} + \phi_2 \mathbb{G}_{ne} \right)^{-1} \label{eq:CT}
\end{align}
$$

We can further simplify Eq. $\eqref{eq:stress_3}$ by defining

$$
\pmb{\varepsilon}^k_\text{rhs} = \pmb{\varepsilon}_{ne}^k - \phi_2 \left( \mathbb{G}_{ne} : \pmb{\sigma}^k + \mathbf{B}_{ne} \right).
\label{eq:eps_rhs}
$$

In this manner, the stress tensor can be expressed as

$$
\pmb{\sigma}^{k+1} = \mathbb{C}_T : \left( \pmb{\varepsilon}_0 + \pmb{\varepsilon}^{k+1} - \pmb{\varepsilon}^k_\text{rhs} \right).
\label{eq:stress_4}
$$

## Linearized momentum balance equation

Finally, the linearized momentum balance equation reads

$$
\nabla \cdot \mathbb{C}_T : \pmb{\varepsilon}^{k+1} = \mathbf{f} + \nabla \cdot \mathbb{C}_T : \left( \pmb{\varepsilon}_\text{rhs}^k - \pmb{\varepsilon}_0 \right).
\label{eq:mom_1}
$$




	
## Weak formulation