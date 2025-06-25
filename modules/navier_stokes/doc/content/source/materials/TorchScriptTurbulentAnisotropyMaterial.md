# TorchScriptTurbulentAnisotropyMaterial

This object computes the Turbulent Anisotropy tensor of a flow using a Neural Network based turbulence closure. For this specific case, the turbulent anisotropy is derived from an Algebraic Reynolds Stress Model from Girimaji [1]. 

## Formulation
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

To make this model mesh with the existing RANS models in MOOSE, the anisotropy tensor must be split into two components, the eddy viscosity and anisotropy correction.
This is because the TKE, TKED, and momentum equations are already eddy viscosity driven. As a result, this material outputs the typical eddy viscosity and an anisotropy correction tensor $b'_{ij}$.

$$
    \mu_t = - \rho \frac{k^2}{\epsilon}G_1 \\
    b'_{ij} = \rho (G_2 (s_{ik}r_{kj} - r_{ik}s_{kj}) + G_3 (s_{ik}s_{kj} - \delta_{ij} s_{nm}s_{nm}))
$$

Thus this material outputs the eddy viscosity through the mu_t_name and the anisotropy corrections with the names property_prefix + "_" + ij. 

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