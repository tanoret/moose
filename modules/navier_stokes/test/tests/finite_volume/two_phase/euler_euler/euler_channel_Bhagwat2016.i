# phase 1 is water, phase 2 is gas
# Bhagwat & Ghajar (2016)
# flow regime slug
# void fraction          0.455    0.55    0.6    0.317    0.417    0.475
# liquid velocity(m/s)   0.15   0.15     0.15    0.45     0.45     0.45
# gas velocity(m/s)      0.21   0.42     0.65    0.39     0.38     0.58

# flow regime bubbly
# void fraction          0.125    0.21    0.27
# mixture velocity(m/s)  1.0501   1.0502  1.0501
mu = 1.002e-3
rho_1 = 998.19
rho_2 = 1.204
pipe_diameter = 0.0127 # m
length = 0.89 # m
U_1 = 0.15   # m/s
U_2 = 0.21   # m/s
g = -9.81
inlet_phase_2 = 0.455
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'

[GlobalParams]
[]

[UserObjects]
  [rc_phase_1]
    type = INSFVRhieChowInterpolator
    u = vel_x_phase_1
    v = vel_y_phase_1
    pressure = pressure
  []
  [rc_phase_2]
    type = INSFVRhieChowInterpolator
    u = vel_x_phase_2
    v = vel_y_phase_2
    pressure = pressure
  []
[]

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = '${length}'
    ymin = '${fparse -pipe_diameter / 2}'
    ymax = '${fparse pipe_diameter / 2}'
    nx = 700
    ny = 10
  []
  uniform_refine = 0
[]

[Variables]
  [vel_x_phase_1]
    type = INSFVVelocityVariable
    initial_condition = ${U_1}
  []
  [vel_y_phase_1]
    type = INSFVVelocityVariable
    initial_condition = 0
  []
  [vel_x_phase_2]
    type = INSFVVelocityVariable
    initial_condition = 1
  []
  [vel_y_phase_2]
    type = INSFVVelocityVariable
    initial_condition = ${U_2}
  []
  [pressure]
    type = INSFVPressureVariable
    two_term_boundary_expansion = false
  []
  [alpha_1_var]
    type = INSFVScalarFieldVariable
    initial_condition = 0.5
  []
  # [pressure_hydro]
  #   type = INSFVPressureVariable
  #   two_term_boundary_expansion = false
  # []
  # [lambda]
  #   family = SCALAR
  #   order = FIRST
  # []
[]

[FVKernels]
  [mass]
    type = WCNSFV2PMassAdvection
    variable = pressure
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho_phases = 'rho_1 rho_2'
    fds = 'alpha_1_var alpha_2_var'
    rhie_chow_user_objects = 'rc_phase_1 rc_phase_2'
  []

  [alpha_1_time]
    type = WCNSFV2PPhaseTimeDerivative
    variable = alpha_1_var
    drho_dt = 0.0
    rho = 'rho_1'
  []
  [alpha_1_advection]
    type = WCNSFV2PPhaseAdvection
    variable = alpha_1_var
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rhie_chow_user_object = 'rc_phase_1'
    rho = 'rho_1'
  []
  [alpha_1_level_set]
    type = WCNSFV2PLevelSet
    variable = alpha_1_var
    gamma = 1.0
    interface_thickness = 0.1
    regularizer_control = 1e-3
  []

  # [mean_alpha]
  #   type = FVIntegralValueConstraint
  #   variable = alpha_1_var
  #   lambda = lambda
  #   phi0 = 0.5
  # []
  # [mean_zero_pressure]
  #   type = FVIntegralValueConstraint
  #   variable = pressure
  #   lambda = lambda
  #   phi0 = 0.0
  # []

  [u_time_1]
    type = WCNSFV2PMomentumTimeDerivative
    variable = vel_x_phase_1
    rhie_chow_user_object = 'rc_phase_1'
    drho_dt = 0.0
    rho = 'rho_1'
    fd = 'alpha_1_var'
    momentum_component = 'x'
  []
  [u_advection_1]
    type = WCNSFV2PMomentumAdvection
    variable = vel_x_phase_1
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rhie_chow_user_object = 'rc_phase_1'
    rho = 'rho_1'
    fd = 'alpha_1_var'
    momentum_component = 'x'
  []
  [u_viscosity_1]
    type = WCNSFV2PMomentumDiffusion
    variable = vel_x_phase_1
    mu = ${mu}
    fd = 'alpha_1_var'
    rhie_chow_user_object = 'rc_phase_1'
    complete_expansion = true
    u = vel_x_phase_1
    v = vel_y_phase_1
    momentum_component = 'x'
  []
  [u_pressure_1]
    type = WCNSFV2PMomentumPressure
    variable = vel_x_phase_1
    fd = 'alpha_1_var'
    rhie_chow_user_object = 'rc_phase_1'
    momentum_component = 'x'
    pressure = pressure
  []
  [u_buoyant_1]
    type = WCNSFV2PMomentumGravity
    variable = vel_x_phase_1
    rho = 'rho_1'
    fd = 'alpha_1_var'
    rho_mixture = 'rho_mixture'
    rhie_chow_user_object = 'rc_phase_1'
    momentum_component = 'x'
    gravity = '${g} 0 0'
  []
  [v_time_1]
    type = WCNSFV2PMomentumTimeDerivative
    variable = vel_y_phase_1
    rhie_chow_user_object = 'rc_phase_1'
    drho_dt = 0.0
    rho = 'rho_1'
    fd = 'alpha_1_var'
    momentum_component = 'y'
  []
  [v_advection_1]
    type = WCNSFV2PMomentumAdvection
    variable = vel_y_phase_1
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rhie_chow_user_object = 'rc_phase_1'
    rho = 'rho_1'
    fd = 'alpha_1_var'
    momentum_component = 'y'
  []
  [v_viscosity_1]
    type = WCNSFV2PMomentumDiffusion
    variable = vel_y_phase_1
    mu = ${mu}
    fd = 'alpha_1_var'
    rhie_chow_user_object = 'rc_phase_1'
    complete_expansion = true
    u = vel_x_phase_1
    v = vel_y_phase_1
    momentum_component = 'y'
  []
  [v_pressure_1]
    type = WCNSFV2PMomentumPressure
    variable = vel_y_phase_1
    fd = 'alpha_1_var'
    rhie_chow_user_object = 'rc_phase_1'
    momentum_component = 'y'
    pressure = pressure
  []
  [v_buoyant_1]
    type = WCNSFV2PMomentumGravity
    variable = vel_y_phase_1
    rho = 'rho_1'
    fd = 'alpha_1_var'
    rho_mixture = 'rho_mixture'
    rhie_chow_user_object = 'rc_phase_1'
    momentum_component = 'y'
    gravity = '${g} 0 0'
  []

  [u_time_2]
    type = WCNSFV2PMomentumTimeDerivative
    variable = vel_x_phase_2
    rhie_chow_user_object = 'rc_phase_2'
    drho_dt = 0.0
    rho = 'rho_2'
    fd = 'alpha_2_var'
    momentum_component = 'x'
  []
  [u_advection_2]
    type = WCNSFV2PMomentumAdvection
    variable = vel_x_phase_2
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rhie_chow_user_object = 'rc_phase_2'
    rho = 'rho_2'
    fd = 'alpha_2_var'
    momentum_component = 'x'
  []
  [u_viscosity_2]
    type = WCNSFV2PMomentumDiffusion
    variable = vel_x_phase_2
    mu = ${mu}
    fd = 'alpha_2_var'
    rhie_chow_user_object = 'rc_phase_2'
    complete_expansion = true
    u = vel_x_phase_2
    v = vel_y_phase_2
    momentum_component = 'x'
  []
  [u_pressure_2]
    type = WCNSFV2PMomentumPressure
    variable = vel_x_phase_2
    fd = 'alpha_2_var'
    rhie_chow_user_object = 'rc_phase_2'
    momentum_component = 'x'
    pressure = pressure
  []
  [u_buoyant_2]
    type = WCNSFV2PMomentumGravity
    variable = vel_x_phase_2
    rho = 'rho_2'
    fd = 'alpha_2_var'
    rho_mixture = 0.0
    rhie_chow_user_object = 'rc_phase_2'
    momentum_component = 'x'
    gravity = '${g} 0 0'
  []
  [v_time_2]
    type = WCNSFV2PMomentumTimeDerivative
    variable = vel_y_phase_2
    rhie_chow_user_object = 'rc_phase_2'
    drho_dt = 0.0
    rho = 'rho_2'
    fd = 'alpha_2_var'
    momentum_component = 'y'
  []
  [v_advection_2]
    type = WCNSFV2PMomentumAdvection
    variable = vel_y_phase_2
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rhie_chow_user_object = 'rc_phase_2'
    rho = 'rho_2'
    fd = 'alpha_2_var'
    momentum_component = 'y'
  []
  [v_viscosity_2]
    type = WCNSFV2PMomentumDiffusion
    variable = vel_y_phase_2
    mu = ${mu}
    fd = 'alpha_2_var'
    rhie_chow_user_object = 'rc_phase_2'
    complete_expansion = true
    u = vel_x_phase_2
    v = vel_y_phase_2
    momentum_component = 'y'
  []
  [v_pressure_2]
    type = WCNSFV2PMomentumPressure
    variable = vel_y_phase_2
    fd = 'alpha_2_var'
    rhie_chow_user_object = 'rc_phase_2'
    momentum_component = 'y'
    pressure = pressure
  []
  [v_buoyant_2]
    type = WCNSFV2PMomentumGravity
    variable = vel_y_phase_2
    rho = 'rho_2'
    fd = 'alpha_2_var'
    rho_mixture = 'rho_mixture'
    rhie_chow_user_object = 'rc_phase_2'
    momentum_component = 'y'
    gravity = '${g} 0 0'
  []

  # [pressure_hydro_laplacian]
  #   type = FVDiffusion
  #   variable = pressure_hydro
  #   coeff = 1.0
  # []
  # [divergence_rho_g_mixture]
  #   type = FVDivergence
  #   variable = pressure_hydro
  #   vector_field = rho_g_mixture
  # []
[]

[FVBCs]
  [inlet-u-1]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_x_phase_1
    functor = '${U_1}'
  []
  [inlet-v-1]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_y_phase_1
    functor = '0'
  []
  [inlet-u-2]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_x_phase_2
    functor = '${U_2}'
  []
  [inlet-v-2]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_y_phase_2
    functor = '0'
  []

  [walls-u-1]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_x_phase_1
    function = 0
  []
  [walls-v-1]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_y_phase_1
    function = 0
  []
  [walls-u-2]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_x_phase_2
    function = 0
  []
  [walls-v-2]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_y_phase_2
    function = 0
  []

  # [outlet_p]
  #   type = INSFVOutletPressureBC
  #   boundary = 'top'
  #   variable = pressure
  #   function = '0.0'
  # []

  # [outlet_p]
  #   type = INSFVAveragePressureValueBC
  #   boundary = 'right'
  #   variable = pressure
  #   lambda = 'lambda'
  # []

  # [fixed_hydro_pressure]
  #   type = FVDirichletBC
  #   boundary = 'top'
  #   variable = pressure_hydro
  #   value = 0.0
  # []

  # [outlet_p_hydro]
  #   type = INSFVOutletPressureBC
  #   boundary = 'right'
  #   variable = pressure
  #   #functor = 'pressure_hydro'
  #   function = '0.0'
  # []

  [outlet_p_hydro]
    type = INSFVOutflowPressureBC
    boundary = 'right'
    variable = pressure
    u = 'vel_x_phase_1'
    v = 'vel_x_phase_2'
    rho = 'rho_mixture'
    mu = ${mu}
  []

  [inlet-alpha_1]
    type = FVDirichletBC
    boundary = 'left'
    variable = alpha_1_var
    value = ${fparse 1-inlet_phase_2}
  []
[]

# [UserObjects]
#   [set_pressure]
#     type = NSPressurePin
#     variable = pressure
#     # pin_type = 'average'
#     # pressure_average = 'pavg'
#     pin_type = 'point-value'
#     point = '0 0 0'
#     execute_on = 'TIMESTEP_END'
#   []
# []

[AuxVariables]
  [alpha_2_var]
    type = MooseVariableFVReal
    initial_condition = 0.5
  []
[]

[AuxKernels]
  [compute_phase_fraction_2]
    type = ParsedAux
    variable = alpha_2_var
    coupled_variables = 'alpha_1_var'
    expression = '1 - alpha_1_var'
  []
[]

[FunctorMaterials]
  [random_properies]
    type = ADGenericFunctorMaterial
    prop_names = 'rho_1     rho_2'
    prop_values = '${rho_1} ${rho_2}'
  []
  [mixing_material]
    type = NSFVMixtureFunctorMaterial
    phase_1_names = 'rho_1'
    phase_2_names = 'rho_2'
    prop_names = 'rho_mixture'
    phase_1_fraction = 'alpha_1_var'
  []
  [rho_g_components]
    type = ADParsedFunctorMaterial
    functor_names = 'rho_mixture'
    property_name = 'rho_g_components'
    expression = 'rho_mixture * (9.81)'
  []
  [rho_g]
    type = ADGenericVectorFunctorMaterial
    prop_names = 'rho_g_mixture'
    prop_values = '0 rho_g_components 0'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_shift_type '
    petsc_options_value = 'lu       NONZERO'
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  # petsc_options_iname = '-pc_type -pc_factor_shift_type -snes_linesearch_dampin'
  # petsc_options_value = 'lu NONZERO                     0.5'
  # dt = 0.1
  # end_time = 10.0
  [TimeSteppers]
    [iter_dt]
      type = IterationAdaptiveDT
      optimal_iterations = 10
      iteration_window = 5
      growth_factor = 2.0
      cutback_factor = 0.25
      dt = 0.0001
    []
    [const_dt]
      type = ConstantDT
      dt = 5
    []
  []
  automatic_scaling = true
  end_time = 1e5
  nl_max_its = 10
  steady_state_detection = true
  steady_state_tolerance = 1e-5
  nl_rel_tol = 1e-03
  nl_abs_tol = 1e-11
  l_max_its = 5
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
   function = '${rho_1} * ${fparse pipe_diameter} * ${U_1}'
    pp_names = ''
  []
  [lin]
    type = NumLinearIterations
  []
  [cum_lin]
    type = CumulativeValuePostprocessor
    postprocessor = lin
  []
  [pavg]
    type = ElementAverageValue
    variable = pressure
  []
  [vg_x]
    type = SideAverageValue
    boundary = 'right'
    variable = 'vel_x_phase_2'
    outputs=none
  []
  [vg_y]
    type = SideAverageValue
    boundary = 'right'
    variable = 'vel_y_phase_2'
    outputs=none
  []
  [vg_value]
    type = ParsedPostprocessor
    expression = 'sqrt(vg_x*vg_x + vg_y*vg_y)*vg_x/abs(vg_x)'
    pp_names = 'vg_x vg_y'
  []
[]

[Outputs]
  print_linear_residuals = true
  print_nonlinear_residuals = true
  [out]
    type = Exodus
    hide = 'Re lin cum_lin'
  []
  [perf]
    type = PerfGraphOutput
  []
[]
