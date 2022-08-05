# This example demonstrates how the mixing length model can be tuned to match an
# established correlation for pressure drop in a smooth circular pipe.

# The primary input parameters for this example are the system Reynolds number
# and the von Karman constant for the mixing length model. These two parameters
# can be changed here:
Re = 1e5
von_karman_const = 0.22

# Note that for this model (using the wall-distance mixing length for the entire
# pipe) different von Karman constants are optimal for different Reynolds
# numbers.

# This model has been non-dimensionalized. The diameter (D), density (rho), and
# bulk velocity (bulk_u) are all considered unity.
D = 1
total_len = ${fparse 20 * D}
rho = 1
bulk_u = 1

# With those parameters set, the viscosity is then computed in order to reach
# the desired Reynolds number.
mu = ${fparse rho * bulk_u * D / Re}

# Here the DeltaP will be evaluted by using a postprocessor to find the pressure
# at a point that is 10 diameters away from the outlet. (The outlet pressure is
# set to zero.)
L = ${fparse 10 * D}

# We will use the McAdams correlation to find the Darcy friction factor. Note
# that this correlation is valid for fully developed flow in smooth circular
# tubes at 3e4 < Re < 1e6.
f = ${fparse 0.316 * Re^(-0.25)}

# The DeltaP can then be computed using this friction factor as,
ref_delta_P = ${fparse f * L / D * rho * bulk_u^2 / 2}

# Regularization constant
alpha = 15.0

# k-Epsilon standard coefficients
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${total_len}
    ymin = 0
    ymax = ${fparse 0.5 * D}
    nx = 200
    ny = 40
    bias_y = ${fparse 1 / 1.2}
  []
  [rename1]
    type = RenameBoundaryGenerator
    input = gen
    old_boundary = 'left'
    new_boundary = 'inlet'
  []
  [rename2]
    type = RenameBoundaryGenerator
    input = rename1
    old_boundary = 'right'
    new_boundary = 'outlet'
  []
  [rename3]
    type = RenameBoundaryGenerator
    input = rename2
    old_boundary = 'bottom'
    new_boundary = 'symmetry'
  []
  [rename4]
    type = RenameBoundaryGenerator
    input = rename3
    old_boundary = 'top'
    new_boundary = 'wall'
  []
[]

[Outputs]
  exodus = true
[]

[Problem]
  kernel_coverage_check = false
  fv_bcs_integrity_check = true
  coord_type = 'RZ'
  rz_coord_axis = 'X'
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  # The upwind and Rhie-Chow interpolation schemes are used here.
  advected_interp_method='upwind'
  velocity_interp_method='rc'
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolator
    u = u_x
    v = u_y
    # a_u = u_x_aux
    # a_v = u_y_aux
    pressure = pressure
  []
[]

[Variables]
  [u_x]
    type = INSFVVelocityVariable
    initial_condition = 1e-6
  []
  [u_y]
    type = INSFVVelocityVariable
    initial_condition = 1e-6
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [k]
    type = INSFVEnergyVariable
    initial_condition = 1.0
  []
  [epsilon]
    type = INSFVEnergyVariable
    initial_condition = 1.0
  []
[]

[AuxVariables]
  [mixing_len]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    rho = ${rho}
  []

  [u_x_advection]
    type = INSFVMomentumAdvection
    variable = u_x
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_x_viscosity]
    type = INSFVMomentumDiffusion
    variable = u_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_x_viscosity_rans]
    type = INSFVMomentumDiffusion
    variable = u_x
    mu = mu_t
    momentum_component = 'x'
  []
  [u_x_pressure]
    type = INSFVMomentumPressure
    variable = u_x
    momentum_component = 'x'
    pressure = pressure
  []

  [u_y_advection]
    type = INSFVMomentumAdvection
    variable = u_y
    rho = ${rho}
    momentum_component = 'y'
  []
  [u_y_viscosity]
    type = INSFVMomentumDiffusion
    variable = u_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [u_y_viscosity_rans]
    type = INSFVMomentumDiffusion
    variable = u_y
    mu = mu_t
    momentum_component = 'y'
  []
  [u_y_pressure]
    type = INSFVMomentumPressure
    variable = u_y
    momentum_component = 'y'
    pressure = pressure
  []

  [TKE_source_sink]
    type = PINSFVTKESourceSink
    variable = k
    u = u_x
    v = u_y
    rho = ${rho}
    mu_t = mu_t
    epsilon = epsilon
    porosity = 1.0
    effective_balance = false
  []
  [TKE_diffusion]
    type = PINSFVTurbulentDiffusion
    variable = k
    mu_t = mu_t
    porosity = 1.0
    turb_coef = ${sigma_k}
    effective_diffusivity = false
  []
  [TKE_diffusion_laminar]
    type = PINSFVTurbulentDiffusion
    variable = k
    mu_t = ${mu}
    porosity = 1.0
    turb_coef = ${sigma_k}
    effective_diffusivity = false
  []
  [TKE_advection]
    type = PINSFVTurbulentAdvection
    variable = k
    velocity_interp_method = 'average'
    advected_interp_method = 'upwind'
    rho = ${rho}
  []

  [TKED_source_sink]
    type = PINSFVTKEDSourceSink
    variable = epsilon
    u = u_x
    v = u_y
    rho = ${rho}
    mu_t = mu_t
    k = k
    porosity = 1.0
    effective_balance = false
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
  []
  [TKED_diffusion]
    type = PINSFVTurbulentDiffusion
    variable = epsilon
    mu_t = mu_t
    porosity = 1.0
    turb_coef = ${sigma_eps}
    effective_diffusivity = false
  []
  [TKED_diffusion_laminar]
    type = PINSFVTurbulentDiffusion
    variable = epsilon
    mu_t = ${mu}
    porosity = 1.0
    turb_coef = ${sigma_eps}
    effective_diffusivity = false
  []
  [TKED_advection]
    type = PINSFVTurbulentAdvection
    variable = epsilon
    velocity_interp_method = 'average'
    advected_interp_method = 'upwind'
    rho = ${rho}
  []
[]

[AuxVariables]
  [mu_t]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1.0
  []
[]

[AuxKernels]
  [compute_mu_t]
    type = kEpsilonViscosity
    variable = mu_t
    k = k
    epsilon = epsilon
    rho = ${rho}
    C_mu = ${C_mu}
  []
[]

[FVBCs]
  [inlet_u_x]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = u_x
    function = ${bulk_u}
  []
  [inlet_u_y]
    type = INSFVInletVelocityBC
    boundary = 'inlet'
    variable = u_y
    function = '0'
  []
  [inlet-k]
    type = FVDirichletBC
    boundary = 'inlet'
    variable = k
    value = 1
  []
  [inlet-eps]
    type = FVDirichletBC
    boundary = 'inlet'
    variable = epsilon
    value = 1
  []

  [walls_u_x]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = u_x
    function = 0
  []
  [walls_u_y]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = u_y
    function = 0
  []
  [walls_k]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = k
    function = 1
  []
  [walls_eps]
    type = INSFVNoSlipWallBC
    boundary = 'wall'
    variable = epsilon
    function = 1
  []

  [sym_u_x]
    type = INSFVSymmetryVelocityBC
    boundary = 'symmetry'
    variable = u_x
    u = u_x
    v = u_y
    mu = ${mu}
    momentum_component = 'x'
  []
  [sym_u_y]
    type = INSFVSymmetryVelocityBC
    boundary = 'symmetry'
    variable = u_y
    u = u_x
    v = u_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [sym_k]
    type = INSFVSymmetryPressureBC
    boundary = 'symmetry'
    variable = k
  []
  [sym_eps]
    type = INSFVSymmetryPressureBC
    boundary = 'symmetry'
    variable = k
  []

  [sym_p]
    type = INSFVSymmetryPressureBC
    boundary = 'symmetry'
    variable = pressure
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'outlet'
    variable = pressure
    function = '0'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type'
  petsc_options_value = 'lu'
  line_search = 'none'
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-12
[]

[Postprocessors]
  [delta_P]
    type = PointValue
    variable = 'pressure'
    point = '${fparse total_len - L} 0 0'
  []
  [reference_delta_P]
    type = Receiver
    default = ${ref_delta_P}
  []
[]
