mu=0.1
rho=1
alpha=1000.0
advected_interp_method='average'
velocity_interp_method='rc'
zero_vel = 0

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
  [u_reg]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1
  []
  [v_reg]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 1e-6
  []
  [pressure]
    type = INSFVPressureVariable
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
  [zero_vel]
    type = PINSFVSuperficialVelocityVariable
    initial_condition = 0.0
  []
[]

[FVKernels]
  inactive = 'mean-pressure'
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u_reg
    v = v_reg
    rho = ${rho}
    porosity = porosity
  []

  [u_advection]
    type = PINSFVMomentumAdvection
    variable = u_reg
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = zero_vel
    v = zero_vel
    rho = ${rho}
    porosity = porosity
    momentum_component = 'x'
  []
  [u_viscosity]
    type = PINSFVMomentumDiffusion
    variable = u_reg
    mu = ${mu}
    porosity = porosity
    momentum_component = 'x'
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = u_reg
    momentum_component = 'x'
    pressure = pressure
    porosity = porosity
  []
  [u_friction]
    type = PINSFVMomentumFriction
    variable = u_reg
    momentum_component = 'x'
    porosity = porosity
    Darcy_name = 'Darcy_coefficient'
    Forchheimer_name = 'Forchheimer_coefficient'
    rho = ${rho}
  []

  # [u_regularization]
  #   type = INSFVMomentumRegularization
  #   variable = u
  #   momentum_component = 'x'
  #   filter_scaling = ${alpha}
  # []
  [u_regularization_source]
    type = INSFVMomentumRegularizationSource
    variable = u
    momentum_component = 'x'
    var_reg = u_reg
  []

  [v_advection]
    type = PINSFVMomentumAdvection
    variable = v_reg
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = zero_vel
    v = zero_vel
    rho = ${rho}
    porosity = porosity
    momentum_component = 'y'
  []
  [v_viscosity]
    type = PINSFVMomentumDiffusion
    variable = v_reg
    mu = ${mu}
    porosity = porosity
    momentum_component = 'y'
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = v_reg
    momentum_component = 'y'
    pressure = pressure
    porosity = porosity
  []
  [v_friction]
    type = PINSFVMomentumFriction
    variable = v_reg
    momentum_component = 'y'
    porosity = porosity
    Darcy_name = 'Darcy_coefficient'
    Forchheimer_name = 'Forchheimer_coefficient'
    rho = ${rho}
  []

  # [v_regularization]
  #   type = INSFVMomentumRegularization
  #   variable = v
  #   momentum_component = 'y'
  #   filter_scaling = ${alpha}
  # []
  [v_regularization_source]
    type = INSFVMomentumRegularizationSource
    variable = v
    momentum_component = 'x'
    var_reg = v_reg
  []

  [mean-pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
    phi0 = 0.01
  []
[]

[FVBCs]
  inactive = 'free-slip-reg-u free-slip-reg-v free-slip-u free-slip-v'
  [inlet-reg-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = u_reg
    function = '1'
  []
  [inlet-reg-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = v_reg
    function = 0
  []
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

  [no-slip-reg-u]
    type = INSFVNoSlipWallBC
    boundary = 'top'
    variable = u_reg
    function = 0
  []
  [no-slip-reg-v]
    type = INSFVNoSlipWallBC
    boundary = 'top'
    variable = v_reg
    function = 0
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
  [free-slip-reg-u]
    type = INSFVNaturalFreeSlipBC
    boundary = 'top'
    variable = u_reg
    momentum_component = 'x'
  []
  [free-slip-reg-v]
    type = INSFVNaturalFreeSlipBC
    boundary = 'top'
    variable = v_reg
    momentum_component = 'y'
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

  [symmetry-reg-u]
    type = PINSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = u_reg
    u = u_reg
    v = v_reg
    mu = ${mu}
    momentum_component = 'x'
  []
  [symmetry-reg-v]
    type = PINSFVSymmetryVelocityBC
    boundary = 'bottom'
    variable = v_reg
    u = u_reg
    v = v_reg
    mu = ${mu}
    momentum_component = 'y'
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
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
    rho = ${rho}
  []
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
    variable = u_reg
    boundary = 'right'
  []
[]

[Outputs]
  exodus = true
[]
