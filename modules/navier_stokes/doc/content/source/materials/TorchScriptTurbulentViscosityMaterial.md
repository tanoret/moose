# TorchScriptTurbulentViscosityMaterial

This object computes the Turbulent viscosity of a flow using a Neural Network based turbulence closure. For this specific case, the turbulent viscosity is derived from an Algebraic Reynolds Stress Model from Girimaji [1]. 

## Formulation

Using the Boussinesq eddy viscosity assumption and incompressibility, one can model the Reynolds stress as:

$$
    - \rho \overline{u_i' u_j'} = \mu_t S_{ij} - \frac{2}{3} \rho k \delta_{ij}
$$

From Girimaji [1],

$$
    s_{ij} = \frac{k}{\epsilon}S_{ij}  \\
    r_{ij} = \frac{k}{\epsilon} R_{ij}  \\
$$

Where 

$$
    S_{ij} = \frac{1}{2} \left(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i} \right) \\
    R_{ij} = \frac{1}{2} \left(\frac{\partial u_i}{\partial x_j} - \frac{\partial u_j}{\partial x_i} \right)
$$

Girimaji focuses on modeling the Reynolds Stress tensor in terms of its tensor basis from Pope [2] using the anisotropy tensor $b_{ij}$.

$$
    b_{ij} = G_1(s_{ij}) + G_2 (s_{ik}r_{kj} - r_{ik}s_{kj}) + G_3 (s_{ik}s_{kj} - \delta_{ij} s_{nm}s_{nm}) \\
    b_{ij} = \frac{\overline{u_i' u_j'}}{2k} - \frac{1}{3} \delta_{ij}
$$

In the ARSM, the values $G_1$ - $G_3$ are a function of scalar invariants $\eta_1$ and $\eta_2$. Defined as:

$$
    \eta_1 = s_{ij} s_{ij} \\
    \eta_2 = r_{ij} r_{ij}
$$

Solving for the Reynolds Stress,

$$
    -\rho\overline{u_i' u_j'} = -2k\rho b_{ij} - \frac{2}{3} k\rho\delta_{ij} 
$$

One can then substitute this result into the Boussinesq assumption and introduce an extra $\frac{k}{\epsilon}$ to yield,

$$
    \mu_t \frac{\epsilon}{k} s_{ij} = -2k\rho[G_1(s_{ij}) + G_2 (s_{ik}r_{kj} - r_{ik}s_{kj}) + G_3 (s_{ik}s_{kj} - \delta_{ij} s_{nm}s_{nm})]
$$

To solve for $\mu_t$ one can then uses a left contraction with $s_{ij}$,

$$
    \mu_t  = -2\frac{k^2}{\epsilon \eta_1}\rho[G_1(\eta_1) + G_2 s_{ij}(s_{ik}r_{kj} - r_{ik}s_{kj}) + G_3 s_{ij}(s_{ik}s_{kj} - \delta_{ij} s_{nm}s_{nm})]
$$

This value is then stored in the property mu_t_torch


## References
<a id="1">[1]</a>
Girimaji, S. (1996)
"Fully Explicit and Self-Consistent Algebraic Reynolds Stress Model"

<a id="2">[2]</a>
Pope, S. B. (1975)
"A more general effective-viscosity hypothesis"



!syntax parameters /Materials/TorchScriptTurbulentViscosityMaterial

!syntax inputs /Materials/TorchScriptTurbulentViscosityMaterial

!syntax children /Materials/TorchScriptTurbulentViscosityMaterial