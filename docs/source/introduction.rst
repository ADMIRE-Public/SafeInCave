Mathematical Formulation
========================

Something here.

Linear momentum balance equation
--------------------------------

To use Lumache, first install it using pip:

.. math::
   :label: eq:mom_0

   \nabla \cdot \pmb{\sigma} = \mathbf{f}

In Eq. :eq:`eq:mom_0`, the stress is calculated as,

.. math::
   :label: eq:stress_0

   \pmb{\sigma} = \mathbb{C}_0 : \pmb{\varepsilon}_{e}

where :math:`\pmb{\varepsilon}_{e}` is the elastic strain tensor and :math:`\mathbb{C}_0` is the linear elastic 4-th order tensor associated to the spring.

The total strain tensor can be represented as

.. math::
   :label: eq:strain_total

   \pmb{\varepsilon} = \pmb{\varepsilon}_{e} + \pmb{\varepsilon}_{ne} = \pmb{\varepsilon}_{e} + \underbrace{\pmb{\varepsilon}_{ve} + \pmb{\varepsilon}_{ie}}_{\pmb{\varepsilon}_{ne}}


.. math::
   :label: eq:stress_1

   \pmb{\sigma} = \mathbb{C}_0^{-1} : \left( \pmb{\varepsilon} - \pmb{\varepsilon}_{ne} \right)

where 

.. math::
   :label: eq:eps_ne

   \pmb{\varepsilon}_{ne} = \sum_{i=1}^{N_{ne}} \pmb{\varepsilon}_{i}

In general, the non-elastic strain rates have a (non-)linear dependency on the stress tensor :math:`\pmb{\sigma}` and, possibly, on internal parameters :math:`\alpha_i`. For example, for an non-elastic element *i*,

.. math::
   :label: eq:eps_ne_sigma_alpha

   \dot{\pmb{\varepsilon}}_{i} = \dot{\pmb{\varepsilon}}_{i} \left( \pmb{\sigma}, \alpha_i \right)



Numerical formulation
=====================

Time integration
----------------

The strain tensor at time :math:`t + \Delta t` of a given non-elastic element :math:`i` can be approximated by

.. math::
   
   \pmb{\varepsilon}_{i}^{t+\Delta t} = \pmb{\varepsilon}^t_{i} + \Delta t \dot{\pmb{\varepsilon}}_{i}^\theta

where :math:`\dot{\pmb{\varepsilon}}_{i}^\theta = \theta \dot{\pmb{\varepsilon}}_{i}^t + (1 - \theta) \dot{\pmb{\varepsilon}}_{i}^{t+\Delta t}`, and :math:`\theta` can be chosen among 0.0, 0.5 and 1.0 for fully implicit, Crank-Nicolson and explicit time integration, respectively. However, the strain rate :math:`\dot{\pmb{\varepsilon}}_{i}^{t+\Delta t}` is unknown and it will be determined in a iterative process. Therefore, we drop the superscript :math:`t+\Delta t` and replace it by :math:`k+1`, where :math:`k` denotes the iterative level. Therefore, the strain of element :math:`i` at iteration :math:`k+1` is

.. math::
   :label: eq:eps_time_integration

   \pmb{\varepsilon}^{k+1}_{i} = \pmb{\varepsilon}^t_{i} + \Delta t \theta \dot{\pmb{\varepsilon}}^t_{i} + \Delta t (1 - \theta) \dot{\pmb{\varepsilon}}^{k+1}_{i}

.. note::

   Keep in mind that both :math:`\pmb{\varepsilon}^t_{i}` and :math:`\dot{\pmb{\varepsilon}}^t_{i}` are known quantities.




Linearized equation
-------------------

Using Taylor expansions to approximate :math:`\dot{\pmb{\varepsilon}}^{k+1}_{i}` in Eq. :eq:`eq:eps_time_integration` it is possible to show that the stress tensor can be computed as

.. math::
   :label: eq:stress_2

   \pmb{\sigma}^{k+1} = \mathbb{C}_T :
    \left[
        \pmb{\varepsilon}^{k+1}
        - \bar{\pmb{\varepsilon}}_{ne}^k
        + \Delta t (1 - \theta)
            \left( 
               \mathbf{B}_{ne}
               + \mathbb{G}_{ne} : \pmb{\sigma}^k
            \right)
    \right]

where

.. math::

   \bar{\pmb{\varepsilon}}_{ne}^k = \sum_{i=1}^{N_{ne}} \left[ \pmb{\varepsilon}_{i}^t + \Delta \theta \dot{\pmb{\varepsilon}}_{i}^t + \Delta t (1 - \theta) \dot{\pmb{\varepsilon}}^{k}_{i} \right],

.. math::

   \mathbb{C}_T = \left[ \mathbb{C}_0^{-1} + \Delta t (1 - \theta) \mathbb{G}_{ne} \right]^{-1},

.. math::

   \mathbb{G}_{ne} = \sum_{i=1}^{N_{ne}} \mathbb{G}_{i} = \sum_{i=1}^{N_{ne}} \left( \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \pmb{\sigma}} - \frac{1}{h_i} \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \alpha_i} \frac{\partial r_i}{\partial \pmb{\sigma}} \right),

.. math::

   \mathbf{B}_{ne} = \sum_{i=1}^{N_{ne}} \mathbf{B}_{i} = \sum_{i=1}^{N_{ne}} \frac{r_i}{h_i} \frac{\partial \dot{\pmb{\varepsilon}}_{i}}{\partial \alpha_i},

.. math::

   \dot{\pmb{\varepsilon}}_{i} = \dot{\pmb{\varepsilon}}_{i}\left( \pmb{\sigma}, \alpha_i \right),

.. math::

   r_{i} = r_{i}\left( \pmb{\sigma}, \alpha_i \right),

.. math::

   h_i = \frac{\partial r_{i}}{\partial \alpha_{i}}.


.. note::

   The residual function :math:`r_{i} = r_{i}\left( \pmb{\sigma}, \alpha_i \right)` is defined based on the evolution equation of :math:`\alpha_i`. Evidently, if the element :math:`i` has no internal parameter, then :math:`r_{i} = 0`.



In this manner, the linearized stress equilibrium equation (Eq. :eq:`eq:mom_0`) can be expressed as

.. math::
   :label: eq:mom_1

   \nabla \cdot \mathbb{C}_T : \pmb{\varepsilon}^{k+1}
    =
    \mathbf{f}
    + \nabla \cdot \mathbb{C}_T : \pmb{\varepsilon}_\text{rhs}^k

where 

.. math::

       \pmb{\varepsilon}_\text{rhs}^k = \bar{\pmb{\varepsilon}}_{ne}^k - \Delta t \left( 1 - \theta \right) \left( \mathbb{G}_{ne} : \pmb{\sigma}^k + \mathbf{B}_{ne} \right)

.. note::

   It is important to note that :math:`\mathbb{G}_{ne}` is a rank-4 tensor, hence the double dot product :math:`:` between :math:`\mathbb{G}_{ne}` and :math:`\pmb{\sigma}^k`. On the other hand, :math:`\mathbf{B}_{ne}` is a rank-2 tensor.


Weak formulation
----------------

Consider a domain :math:`\Omega` bounded by a surface :math:`\Gamma` outward oriented by a normal vector :math:`\mathbf{n}`. Additionally, consider a vector test function :math:`\mathbf{v} \in \mathcal{V}`, where :math:`\mathcal{V}` is the test function space generated by continuous piecewise linear polynomials. In this manner, the weak formulation of the linearized momentum balance equation can be expressed as, 

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

where :math:`a\left( \mathbf{u}, \mathbf{v} \right)` and :math:`L\left( \mathbf{v} \right)` represent the well-known bilinear a linear operators, respectively.


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

.. code-block:: none

   BEGIN
      INPUTS param1, param2
      IF param1 > param2 THEN
         :math:`\sqrt{param1^2 + param2^2}`
         RETURN param1 - param2
      ELSE
         RETURN param2 - param1
      ENDIF
   END

