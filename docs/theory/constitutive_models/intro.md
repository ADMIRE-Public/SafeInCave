# Introduction

In general, a constitutive model can be represented as illustrated in [](#fig-constitutive-model-1), which shows a serial arrangement of different types of elements (springs, dashpots, etc). The total strain $\pmb{\varepsilon}$ is given by the sum of the individual strain of all elements composing the constitutive model. In this text, we make a distinction between **elastic**, **non-elastic**, and **thermal** strains.

![General representation of a constitutive model.](../../images/constitutive_model_1.png){#fig-constitutive-model-1 width="80%"}

Elastic strains $\pmb{\varepsilon}_{e}$ refer exclusively to time-independent (instantaneous) elastic strains -- in other words, it only includes the deformation of the yellow spring in [](#fig-constitutive-model-1).

Non-elastic strains comprise the viscoelastic ($\pmb{\varepsilon}_{ve}$) and inelastic ($\pmb{\varepsilon}_{ie}$) strains. In the SafeInCave simulator, the only viscoelastic element implemented is the Kelvin-Voigt element, which is described below. For inelastic elements, the SafeInCave simulator provides three options: a viscoplastic element, a dislocation creep element, and a pressure solution creep element.

All these different types of elements can be arbitrarily combined as illustrated in [](#fig-constitutive-model-1).

From the discussion above and from [](#fig-constitutive-model-1), it follows that the total strain can be written as

$$
\pmb{\varepsilon} = \pmb{\varepsilon}_e + \underbrace{\pmb{\varepsilon}_{ve} + \pmb{\varepsilon}_{ie}}_{\pmb{\varepsilon}_{ne}} + \pmb{\varepsilon}_{th}
\quad \rightarrow \quad
\pmb{\varepsilon} = \pmb{\varepsilon}_e + \pmb{\varepsilon}_{ne} + \pmb{\varepsilon}_{th}.
$$

The constitutive models available in the literature, some of which are presented in the following sections, differ from each other by how they consider the **non-elastic** strain tensor $\pmb{\varepsilon}_{ne}$.