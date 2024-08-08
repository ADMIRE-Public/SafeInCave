Mathematical Formulation
========================

Something here.

Linear momentum balance equation
--------------------------------

The linear momentum balance equation considering quasi-static loads can be written as

.. math::
   :label: eq:mom_0

   \nabla \cdot \pmb{\sigma} = \mathbf{f}

with :math:`\mathbf{f}` representing the body forces. In Eq. :eq:`eq:mom_0`, the stress is calculated as,

.. math::
   :label: eq:stress_0

   \pmb{\sigma} = \mathbb{C}_0 : \pmb{\varepsilon}_{e}

where :math:`\pmb{\varepsilon}_{e}` is the elastic strain tensor and :math:`\mathbb{C}_0` is the 4th-order tensor associated to the linear elastic response of the material. In most constitutive models for geomaterials, non-elastic deformations are also present.

.. note::
   In the present work, non-elastic deformation includes all types of deformation that are not instantaneously elastic, that is, viscoelastic (time dependent elastic) and inelastic (viscoplastic, plastic, creep, etc) deformations.

The total strain tensor can be represented as

.. math::
   :label: eq:strain_total

   \pmb{\varepsilon} = \pmb{\varepsilon}_{e} + \pmb{\varepsilon}_{ne} = \pmb{\varepsilon}_{e} + \underbrace{\pmb{\varepsilon}_{ve} + \pmb{\varepsilon}_{ie}}_{\pmb{\varepsilon}_{ne}}

where :math:`\pmb{\varepsilon}_{ve}` and :math:`\pmb{\varepsilon}_{ie}` are the viscoelastic and inelastic strains, respectively, and

.. math::
   :label: eq:eps_ne

   \pmb{\varepsilon}_{ne} = \sum_{i=1}^{N_{ne}} \pmb{\varepsilon}_{i}

with :math:`N_{ne}` denoting the number of non-elastic elements included in the constitutive model. In this manner, the stress tensor can be expressed as

.. math::
   :label: eq:stress_1

   \pmb{\sigma} = \mathbb{C}_0^{-1} : \left( \pmb{\varepsilon} - \pmb{\varepsilon}_{ne} \right)

In general, the non-elastic strain rates have a (non-)linear dependency on the stress tensor :math:`\pmb{\sigma}` and, possibly, on internal parameters :math:`\alpha_i`. For example, for an non-elastic element *i*,

.. math::
   :label: eq:eps_ne_sigma_alpha

   \dot{\pmb{\varepsilon}}_{i} = \dot{\pmb{\varepsilon}}_{i} \left( \pmb{\sigma}, \alpha_i \right)

The circular dependency of the non-elastic strains on the stress tensor :math:`\pmb{\sigma}` makes of Eq. :eq:`eq:mom_0` a non-linear equation. The numerical procedure for treating this non-linearity and solving Eq. :eq:`eq:mom_0` is described below.



Numerical formulation
=====================

Time integration
----------------

The strain tensor at time :math:`t + \Delta t` of a given non-elastic element :math:`i` can be approximated by

.. math::
   
   \pmb{\varepsilon}_{i}^{t+\Delta t} = \pmb{\varepsilon}^t_{i} + \Delta t \dot{\pmb{\varepsilon}}_{i}^\theta

where :math:`\dot{\pmb{\varepsilon}}_{i}^\theta = \theta \dot{\pmb{\varepsilon}}_{i}^t + (1 - \theta) \dot{\pmb{\varepsilon}}_{i}^{t+\Delta t}`, and :math:`\theta` can be chosen among 0.0, 0.5 and 1.0 for fully implicit, Crank-Nicolson and explicit time integration, respectively. However, the strain rate :math:`\dot{\pmb{\varepsilon}}_{i}^{t+\Delta t}` is unknown and it will be determined in a iterative process, so we drop the superscript :math:`t+\Delta t` and replace it by :math:`k+1`, where :math:`k` denotes the iterative level. In this manner, the strain of element :math:`i` at iteration :math:`k+1` is

.. math::
   :label: eq:eps_time_integration

   \pmb{\varepsilon}^{k+1}_{i} = \pmb{\varepsilon}^t_{i} + \Delta t \theta \dot{\pmb{\varepsilon}}^t_{i} + \Delta t (1 - \theta) \dot{\pmb{\varepsilon}}^{k+1}_{i}.

For conciseness, let us consider :math:`\phi_1 = \Delta t \theta` and :math:`\phi_2 = \Delta t (1 - \theta)`. Recalling Eq. :eq:`eq:eps_ne` and substituting Eq. :eq:`eq:eps_time_integration` into Eq. :eq:`eq:stress_1`, the stress tensor at iteration :math:`k+1` is expressed as

.. math::
   :label: eq:stress_2
   
   \pmb{\sigma}^{k+1} = \mathbb{C}_0 : \left( \pmb{\varepsilon}^{k+1} - \pmb{\varepsilon}^t_{ne} - \phi_1 \dot{\pmb{\varepsilon}}^t_{ne} - \phi_2 \dot{\pmb{\varepsilon}}^{k+1}_{ne} \right)

where :math:`\dot{\pmb{\varepsilon}}^{k+1}_{ne}` is obviously unknown, which requires a linearization method for its evaluation.

.. note::

   Keep in mind that both :math:`\pmb{\varepsilon}^t_{i}` and :math:`\dot{\pmb{\varepsilon}}^t_{i}` are known quantities.


Picard's method
---------------

One alternative to linearize Eq. :eq:`eq:stress_2` is to simply consider

.. math::
   
   \dot{\pmb{\varepsilon}}^{k+1}_{i} = \dot{\pmb{\varepsilon}}^{k}_{i}

in which :math:`\dot{\pmb{\varepsilon}}^{k}_{i} = \dot{\pmb{\varepsilon}}_{i} \left( \pmb{\sigma}^k, \alpha^k_i \right)`. As a consequence, the stress tensor is linearized as

.. math::
   
   \pmb{\sigma}^{k+1} = \mathbb{C}_0 : \left( \pmb{\varepsilon}^{k+1} - \pmb{\varepsilon}^t_{ne} - \phi_1 \dot{\pmb{\varepsilon}}^t_{ne} - \phi_2 \dot{\pmb{\varepsilon}}^{k}_{ne} \right)

and the momentum balance equation becomes

.. math::
   :label: eq:mom_picard

   \nabla \cdot \mathbb{C}_0 : \pmb{\varepsilon}^{k+1} = \mathbf{f} + \nabla \cdot \mathbb{C}_0 : \left( \pmb{\varepsilon}^t_{ne} + \phi_1 \dot{\pmb{\varepsilon}}^t_{ne} + \phi_2 \dot{\pmb{\varepsilon}}^{k}_{ne} \right).

Although very simple, Eq. :eq:`eq:mom_picard` requires many iterations to converge and it is often unstable, especially when highly non-linear deformations are present.


Newton's method
---------------

Alternatively, the strain rate can be expanded from iteration :math:`k` to iteration :math:`k+1` by using Taylor series, that is,

.. math::
   :label: eq:eps_newton_0

   \dot{\pmb{\varepsilon}}^{k+1}_{i} = \dot{\pmb{\varepsilon}}^{k}_{i} + \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \pmb{\sigma}} : \delta \pmb{\sigma} + \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \alpha_i} \delta \alpha_i

where :math:`\delta \pmb{\sigma} = \pmb{\sigma}^{k+1} - \pmb{\sigma}^k` and :math:`\delta \alpha_i = \alpha_i^{k+1} - \alpha_i^k`.

.. note::

   The term :math:`\frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \pmb{\sigma}}` is a rank-4 tensor, whereas :math:`\delta \pmb{\sigma}` is a rank-2 tensor, hence the double dot product between them, which results a rank-2 tensor.

The increment of internal variable :math:`\delta \alpha_i` can be obtained by defining a residual equation of the evolution equation of :math:`\alpha_i` and using Newton-Raphson to drive the residue to zero. Considering the residual equation is of the form :math:`r_i = r_i(\pmb{\sigma}, \alpha_i)`, it follows that

.. math::
   :label: eq:delta_alpha

   r_i^{k+1} = r^k + \frac{\partial r_i}{\partial \pmb{\sigma}} : \delta \pmb{\sigma} + \underbrace{\frac{\partial r_i}{\partial \alpha_i}}_{h_i} \delta \alpha_i = 0
   \quad \rightarrow \quad
   \delta \alpha_i = - \frac{1}{h_i} \left( r_i^k + \frac{\partial r_i}{\partial \pmb{\sigma}} : \delta \pmb{\sigma} \right).

Substituting Eq. :eq:`eq:delta_alpha` into Eq. :eq:`eq:eps_newton_0` yields

.. math::

   \dot{\pmb{\varepsilon}}^{k+1}_{i} 
   = \dot{\pmb{\varepsilon}}^{k}_{i} 
   + \underbrace{
      \left( \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \pmb{\sigma}} - \frac{1}{h_i} \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \alpha_i} \frac{\partial r_i}{\partial \pmb{\sigma}} \right)
      }_{\mathbb{G}_i} : \delta \pmb{\sigma} 
   - \underbrace{ \frac{r_i^k}{h_i} \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \alpha_i} }_{\mathbf{B}_i}

Considering all non-elastic elements,   

.. math::
   :label: eq:eps_newton_1

   \dot{\pmb{\varepsilon}}^{k+1}_{ne} 
   = \dot{\pmb{\varepsilon}}^{k}_{ne}
   + \mathbb{G}_{ne} : \delta \pmb{\sigma} 
   - \mathbf{B}_{ne}

where :math:`\mathbb{G}_{ne} = \sum_{i=1}^{N_{ne}} \mathbb{G}_i` and :math:`\mathbf{B}_{ne} = \sum_{i=1}^{N_{ne}} \mathbf{B}_i`.

Finally, substituting Eq. :eq:`eq:eps_newton_1` into Eq. :eq:`eq:stress_2` leads to

.. math::
   :label: eq:stress_3

   \pmb{\sigma}^{k+1} = \mathbb{C}_T : \left[
      \pmb{\varepsilon}^{k+1}
      - \bar{\pmb{\varepsilon}}^k_{ne}
      + \phi_2 \left(
         \mathbb{G}_{ne} : \pmb{\sigma}^k
         + \mathbf{B}_{ne}
      \right)
   \right]

where :math:`\bar{\pmb{\varepsilon}}^k_{ne} = \pmb{\varepsilon}^t_{ne} + \phi_1 \dot{\pmb{\varepsilon}}^t_{ne} + \phi_2 \dot{\pmb{\varepsilon}}^{k}_{ne}` and the consistent tangent matrix :math:`\mathbb{C}_T` is given by

.. math::

   \mathbb{C}_T = \left( \mathbb{C}_0^{-1} + \phi_1 \mathbb{G}_{ne} \right)^{-1}.

We can further simplify Eq. :eq:`eq:stress_3` by defining

.. math::

   \pmb{\varepsilon}_{\text{rhs}}^k = \bar{\pmb{\varepsilon}}^k_{ne} - \phi_2 \left(
         \mathbb{G}_{ne} : \pmb{\sigma}^k
         + \mathbf{B}_{ne} \right)

In this manner, the stress tensor can be expressed as

.. math::
   :label: eq:stress_4

   \pmb{\sigma}^{k+1} = \mathbb{C}_T : \left(
      \pmb{\varepsilon}^{k+1}
      - \pmb{\varepsilon}^k_{\text{rhs}}
   \right)

and the linearized momentum balance equation becomes

.. math::
   :label: eq:mom_2

   \nabla \cdot \mathbb{C}_T : \pmb{\varepsilon}^{k+1}
    =
    \mathbf{f}
    + \nabla \cdot \mathbb{C}_T : \pmb{\varepsilon}_\text{rhs}^k

   


.. note::

   It is important to note that :math:`\mathbb{G}_{ne}` is a rank-4 tensor, hence the double dot product :math:`:` between :math:`\mathbb{G}_{ne}` and :math:`\pmb{\sigma}^k`. On the other hand, :math:`\mathbf{B}_{ne}` is a rank-2 tensor.


Weak formulation
----------------

Consider a domain :math:`\Omega` bounded by a surface :math:`\Gamma` outward oriented by a normal vector :math:`\mathbf{n}`. Additionally, consider a vector **test** function :math:`\mathbf{v} \in \mathcal{V}` and a vector **trial** function :math:`\mathbf{u} \in \mathcal{V}`, where :math:`\mathcal{V}` is the test function space generated by continuous piecewise linear polynomials. In this manner, the weak formulation of the linearized momentum balance equation can be expressed as, 

.. math::

   \underbrace{
        \int_\Omega \mathbb{C}_T : \pmb{\varepsilon} \left( \mathbf{u}^{k+1} \right) : \pmb{\varepsilon} \left( \mathbf{v} \right) \text{d} \Omega
    }_{
        a\left( \mathbf{u}, \mathbf{v} \right)
    }
    =
    \underbrace{
        \int_\Omega \mathbf{f} \cdot \mathbf{v} \text{d} \Omega
        +
        \int_\Gamma \mathbf{t} \cdot \mathbf{v} \text{d} \Gamma
        +
        \int_\Omega \mathbb{C}_T : \pmb{\varepsilon}_\text{rhs}^k : \pmb{\varepsilon} \left( \mathbf{v} \right) \text{d} \Omega
    }_{
        L\left( \mathbf{v} \right)
    }

where :math:`a\left( \mathbf{u}, \mathbf{v} \right)` and :math:`L\left( \mathbf{v} \right)` represent the well-known bilinear a linear operators, respectively. The term :math:`\mathbf{t}` is the traction vector applied at the portion of :math:`\Gamma` where Neumann boundary conditions are imposed. Additionally, due to small strain assumption,

.. math::

   \pmb{\varepsilon}(\mathbf{w}) = \frac{1}{2} \left( \nabla \mathbf{w} + \nabla \mathbf{w}^T \right),

in which :math:`\mathbf{w} \in \mathcal{V}`.



.. _constitutive-models-section:

Constitutive models
===================


Viscoelastic element
--------------------

.. math::
   :label: eq:eps_rate_ve_0

   \pmb{\sigma} = \underbrace{\mathbb{C}_1 : \pmb{\varepsilon}_{ve}}_{\text{spring}} + \underbrace{\eta_1 \dot{\pmb{\varepsilon}}_{ve}}_{\text{dashpot}}
    \quad \Rightarrow \quad
    \dot{\pmb{\varepsilon}}_{ve} = \frac{1}{\eta_1} \left( \pmb{\sigma} - \mathbb{C}_1 : \pmb{\varepsilon}_{ve} \right)

Dislocation creep element
-------------------------

.. math::
   :label: eq:eps_rate_dc_0

   \dot{\pmb{\varepsilon}}_{cr} = A \exp \left( -\frac{Q}{RT} \right) q^{n-1} \mathbf{s}

Viscoplastic element
--------------------

.. math::
   :label: eq:eps_rate_vp_0

   \dot{\pmb{\varepsilon}}_{vp} = \mu_1 \left\langle \dfrac{ F_{vp} }{F_0} \right\rangle^{N_1} \dfrac{\partial F_{vp}}{\partial \pmb{\sigma}}

.. math::
   :label: eq:F_vp_0

   F_{vp}(\pmb{\sigma}, \alpha) = J_2 - (-\alpha I_1^{n} + \gamma I_1^2) \left[ \exp{(\beta_1 I_1)} - \beta \cos(3\theta) \right]^m

.. math::
   :label: eq:alpha_0

   \alpha = a_1 \left[ \left( \frac{a_1}{\alpha_0} \right)^{1/\eta} + \xi \right]^{-\eta}, \quad \text{where} \quad \xi = \int_{t_0}^t \sqrt{ \dot{\pmb{\varepsilon}}_{vp} : \dot{\pmb{\varepsilon}}_{vp} } \mathrm{dt}

Algorithms
~~~~~~~~~~

.. code-block:: latex

   BEGIN
      INPUTS param1, param2
      IF param1 > param2 THEN
         :math:`\sqrt{param1^2 + param2^2}`
         RETURN param1 - param2
      ELSE
         RETURN param2 - param1
      ENDIF
   END


