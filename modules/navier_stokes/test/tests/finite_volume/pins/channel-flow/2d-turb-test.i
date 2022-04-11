mu=1e-5
rho=1
advected_interp_method='average'
velocity_interp_method='rc'

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
    xmax = 5
    ymin = 0
    ymax = 1
    nx = 40
    ny = 20
  []
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[UserObjects]
  [rc]
    type = PINSFVRhieChowInterpolator
    u = u
    v = v
    pressure = pressure
    porosity = porosity
  []
[]

[Variables]
  inactive = 'lambda'
  [u]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1
  []
  [v]
    type = PINSFVSuperficialVelocityVariable
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
  [lambda]
    family = SCALAR
    order = FIRST
  []
[]

[AuxVariables]
  [porosity]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 0.5
  []
  [mu_t]
    family = MONOMIAL
    order = CONSTANT
    fv = true
    initial_condition = 1.0
  []
[]

[FVKernels]
  inactive = 'mean-pressure u_friction v_friction'
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
  []

  [u_advection]
    type = PINSFVMomentumAdvection
    variable = u
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    porosity = porosity
    momentum_component = 'x'
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = u
    mu = ${mu}
    porosity = porosity
    momentum_component = 'x'
  []
  [u_viscosity_turbulent]
    type = PINSFVMomentumDiffusion
    variable = u
    mu = mu_t
    porosity = porosity
    momentum_component = 'x'
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
    porosity = porosity
  []
  [u_friction]
    type = PINSFVMomentumFriction
    variable = u
    momentum_component = 'x'
    porosity = porosity
    Darcy_name = 'Darcy_coefficient'
    Forchheimer_name = 'Forchheimer_coefficient'
    rho = ${rho}
  []

  [v_advection]
    type = PINSFVMomentumAdvection
    variable = v
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    porosity = porosity
    momentum_component = 'y'
  []
  [v_viscosity]
    type = PINSFVMomentumDiffusion
    variable = v
    mu = ${mu}
    porosity = porosity
    momentum_component = 'y'
  []
  [v_viscosity_turbulent]
    type = PINSFVMomentumDiffusion
    variable = v
    mu = mu_t
    porosity = porosity
    momentum_component = 'y'
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
    porosity = porosity
  []
  [v_friction]
    type = PINSFVMomentumFriction
    variable = v
    momentum_component = 'y'
    porosity = porosity
    Darcy_name = 'Darcy_coefficient'
    Forchheimer_name = 'Forchheimer_coefficient'
    rho = ${rho}
  []

  [TKE_source_sink]
    type = PINSFVTKESourceSink
    variable = k
    u = u
    v = v
    rho = ${rho}
    mu_t = mu_t
    epsilon = epsilon
    porosity = porosity
    effective_balance = false
  []
  [TKE_diffusion]
    type = PINSFVTurbulentDiffusion
    variable = k
    mu_t = mu_t
    porosity = porosity
    turb_coef = ${sigma_k}
    effective_diffusivity = false
  []
  [TKE_diffusion_laminar]
    type = PINSFVTurbulentDiffusion
    variable = k
    mu_t = ${mu}
    porosity = porosity
    turb_coef = ${sigma_k}
    effective_diffusivity = false
  []
  [TKE_advection]
    type = PINSFVTurbulentAdvection
    variable = k
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
  []

  [TKED_source_sink]
    type = PINSFVTKEDSourceSink
    variable = epsilon
    u = u
    v = v
    rho = ${rho}
    mu_t = mu_t
    k = k
    porosity = porosity
    effective_balance = false
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
  []
  [TKED_diffusion]
    type = PINSFVTurbulentDiffusion
    variable = epsilon
    mu_t = mu_t
    porosity = porosity
    turb_coef = ${sigma_eps}
    effective_diffusivity = false
  []
  [TKED_diffusion_laminar]
    type = PINSFVTurbulentDiffusion
    variable = epsilon
    mu_t = ${mu}
    porosity = porosity
    turb_coef = ${sigma_eps}
    effective_diffusivity = false
  []
  [TKED_advection]
    type = PINSFVTurbulentAdvection
    variable = epsilon
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    rho = ${rho}
  []

  [mean-pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
    phi0 = 0.01
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
  inactive = 'free-slip-u free-slip-v free-slip-k free-slip-eps'
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = u
    function = '1'
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = v
    function = 0
  []
  [inlet-k]
    type = FVDirichletBC
    boundary = 'left'
    variable = k
    value = 1
  []
  [inlet-eps]
    type = FVDirichletBC
    boundary = 'left'
    variable = epsilon
    value = 1
  []

  [no-slip-u]
    type = INSFVNoSlipWallBC
    boundary = 'top'
    variable = u
    function = 0
  []
  [no-slip-v]
    type = INSFVNoSlipWallBC
    boundary = 'top'
    variable = v
    function = 0
  []
  [no-slip-k]
    type = INSFVNoSlipWallBC
    boundary = 'top'
    variable = k
    function = '1'
  []
  [no-slip-eps]
    type = INSFVNoSlipWallBC
    boundary = 'top'
    variable = epsilon
    function = '1'
  []

  [free-slip-u]
    type = INSFVNaturalFreeSlipBC
    boundary = 'top'
    variable = u
    momentum_component = 'x'
  []
  [free-slip-v]
    type = INSFVNaturalFreeSlipBC
    boundary = 'top'
    variable = v
    momentum_component = 'y'
  []
  [free-slip-k]
    type = INSFVNaturalFreeSlipBC
    boundary = 'top'
    variable = k
    function = '1'
  []
  [free-slip-eps]
    type = INSFVNaturalFreeSlipBC
    boundary = 'top'
    variable = epsilon
    function = '1'
  []

  [symmetry-u]
    type = PINSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = u
    u = u
    v = v
    mu = ${mu}
    momentum_component = 'x'
  []
  [symmetry-v]
    type = PINSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = v
    u = u
    v = v
    mu = ${mu}
    momentum_component = 'y'
  []
  [symmetry-p]
    type = INSFVSymmetryPressureBC
    boundary = 'bottom'
    variable = pressure
  []

  [outlet-p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = 0
  []
[]

[Materials]
  [darcy]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'Darcy_coefficient Forchheimer_coefficient'
    prop_values = '0.1 0.1 0.1 0.1 0.1 0.1'
  []
[]

[Executioner]
  type = Steady
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm      200                lu           NONZERO'
  line_search = 'none'
  nl_rel_tol = 1e-11
  nl_abs_tol = 1e-14
[]

# Some basic Postprocessors to visually examine the solution
[Postprocessors]
  [inlet-p]
    type = SideAverageValue
    variable = pressure
    boundary = 'left'
  []
  [outlet-u]
    type = SideIntegralVariablePostprocessor
    variable = u
    boundary = 'right'
  []
[]

[Outputs]
  exodus = true
[]
