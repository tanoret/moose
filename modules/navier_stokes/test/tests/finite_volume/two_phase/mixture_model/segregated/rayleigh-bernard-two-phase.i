mu = 1.0
rho = 1e3
mu_d = 0.3
rho_d = 1.0
dp = 0.01
U_lid = 0.0
g = -9.81
advected_interp_method = 'upwind'

# Currently required
k = 1
k_d = 1
cp = 1
cp_d = 1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = .1
    ymin = 0
    ymax = .1
    nx = 11
    ny = 11
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system phi_system'
[]

[AuxVariables]
  [U]
    type = MooseVariableFVReal
    family = MONOMIAL
    fv = true
    order = CONSTANT
  []
[]

[Executioner/TimeIntegrators]
  [ImplicitEuler]
    type = ImplicitEuler
  []
[]

[AuxKernels]
  [mixture_phase_1_fractionphase_1_out]
    type = FunctorMaterialRealAux
    block = ANY_BLOCK_ID
    execute_on = 'INITIAL TIMESTEP_END'
    functor = phase_1
    variable = phase_1_out
  []
  [mixture_mixture_materialcp_mixture_out]
    type = FunctorMaterialRealAux
    block = ANY_BLOCK_ID
    execute_on = 'INITIAL TIMESTEP_END'
    functor = cp_mixture
    variable = cp_mixture_out
  []
  [mixture_mixture_materialk_mixture_out]
    type = FunctorMaterialRealAux
    block = ANY_BLOCK_ID
    execute_on = 'INITIAL TIMESTEP_END'
    functor = k_mixture
    variable = k_mixture_out
  []
  [mixture_mixture_materialmu_mixture_out]
    type = FunctorMaterialRealAux
    block = ANY_BLOCK_ID
    execute_on = 'INITIAL TIMESTEP_END'
    functor = mu_mixture
    variable = mu_mixture_out
  []
  [mixture_mixture_materialrho_mixture_out]
    type = FunctorMaterialRealAux
    block = ANY_BLOCK_ID
    execute_on = 'INITIAL TIMESTEP_END'
    functor = rho_mixture
    variable = rho_mixture_out
  []
  [mixture_dispersed_dragDarcy_coefficient_out]
    type = FunctorMaterialRealAux
    block = ANY_BLOCK_ID
    execute_on = 'INITIAL TIMESTEP_END'
    functor = Darcy_coefficient
    variable = Darcy_coefficient_out
  []
  [mixture_dispersed_dragDarcy_coefficient_vec_out_x]
    type = FunctorMaterialRealVectorValueAux
    block = ANY_BLOCK_ID
    component = 0
    execute_on = 'INITIAL TIMESTEP_END'
    functor = Darcy_coefficient_vec
    variable = Darcy_coefficient_vec_out_x
  []
  [mixture_dispersed_dragDarcy_coefficient_vec_out_y]
    type = FunctorMaterialRealVectorValueAux
    block = ANY_BLOCK_ID
    component = 1
    execute_on = 'INITIAL TIMESTEP_END'
    functor = Darcy_coefficient_vec
    variable = Darcy_coefficient_vec_out_y
  []
  [mixture_dispersed_dragDarcy_coefficient_vec_out_z]
    type = FunctorMaterialRealVectorValueAux
    block = ANY_BLOCK_ID
    component = 2
    execute_on = 'INITIAL TIMESTEP_END'
    functor = Darcy_coefficient_vec
    variable = Darcy_coefficient_vec_out_z
  []
[]

[AuxVariables]
  [Darcy_coefficient_out]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [Darcy_coefficient_vec_out_x]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [Darcy_coefficient_vec_out_y]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [Darcy_coefficient_vec_out_z]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [cp_mixture_out]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [k_mixture_out]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [mu_mixture_out]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [phase_1_out]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
  [rho_mixture_out]
    type = MooseVariableConstMonomial
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Functions]
  [1e-12]
    type = ConstantFunction
    value = 1e-12
  []
  [0.2]
    type = ConstantFunction
    value = 0.2
  []
[]

[ICs]
  [flow_vel_x_ic]
    type = FunctionIC
    function = 1e-12
    variable = vel_x
  []
  [flow_vel_y_ic]
    type = FunctionIC
    function = 1e-12
    variable = vel_y
  []
  [flow_pressure_ic]
    type = FunctionIC
    function = 0.2
    variable = pressure
  []
[]

[LinearFVBCs]
  [vel_x_bottom]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    functor = 0
    variable = vel_x
  []
  [vel_y_bottom]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = bottom
    functor = 0
    variable = vel_y
  []
  [vel_x_left]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = left
    functor = 0
    variable = vel_x
  []
  [vel_y_left]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = left
    functor = 0
    variable = vel_y
  []
  [vel_x_right]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = right
    functor = 0
    variable = vel_x
  []
  [vel_y_right]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = right
    functor = 0
    variable = vel_y
  []
  [vel_x_top]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = top
    functor = 0.0
    variable = vel_x
  []
  [vel_y_top]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = top
    functor = 0
    variable = vel_y
  []
  [pressure_extrapolation_top_left_right_bottom]
    type = LinearFVExtrapolatedPressureBC
    boundary = 'top left right bottom'
    use_two_term_expansion = true
    variable = pressure
  []
[]

[LinearFVKernels]
  [flow_p_diffusion]
    type = LinearFVAnisotropicDiffusion
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
    variable = pressure
  []
  [flow_HbyA_divergence]
    type = LinearFVDivergence
    face_flux = HbyA
    force_boundary_execution = true
    variable = pressure
  []
  [flow_ins_momentum_time_x]
    type = LinearFVTimeDerivative
    factor = rho_mixture
    variable = vel_x
  []
  [flow_ins_momentum_time_y]
    type = LinearFVTimeDerivative
    factor = rho_mixture
    variable = vel_y
  []
  [flow_ins_momentum_flux_x]
    type = LinearWCNSFVMomentumFlux
    advected_interp_method = upwind
    momentum_component = x
    mu = mu_mixture
    rhie_chow_user_object = ins_rhie_chow_interpolator
    u = vel_x
    use_deviatoric_terms = false
    use_nonorthogonal_correction = false
    v = vel_y
    variable = vel_x
  []
  [flow_ins_momentum_flux_y]
    type = LinearWCNSFVMomentumFlux
    advected_interp_method = upwind
    momentum_component = y
    mu = mu_mixture
    rhie_chow_user_object = ins_rhie_chow_interpolator
    u = vel_x
    use_deviatoric_terms = false
    use_nonorthogonal_correction = false
    v = vel_y
    variable = vel_y
  []
  [flow_ins_momentum_pressure_x]
    type = LinearFVMomentumPressure
    momentum_component = x
    pressure = pressure
    variable = vel_x
  []
  [flow_ins_momentum_pressure_y]
    type = LinearFVMomentumPressure
    momentum_component = y
    pressure = pressure
    variable = vel_y
  []
  [flow_ins_momentum_gravity_y]
    type = LinearFVSource
    source_density = -9.810000
    variable = vel_y
  []
[]

[Materials]
  [flow_ins_speed_material]
    type = ADVectorMagnitudeFunctorMaterial
    execute_on = ALWAYS
    outputs = none
    vector_magnitude_name = speed
    x_functor = vel_x
    y_functor = vel_y
  []
[]

[UserObjects]
  [ins_rhie_chow_interpolator]
    type = RhieChowMassFlux
    p_diffusion_kernel = flow_p_diffusion
    pressure = pressure
    rho = rho_mixture
    u = vel_x
    v = vel_y
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    family = MONOMIAL
    fv = true
    order = CONSTANT
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    family = MONOMIAL
    fv = true
    order = CONSTANT
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    family = MONOMIAL
    fv = true
    order = CONSTANT
    solver_sys = pressure_system
  []
[]

[Functions]
  [0]
    type = ConstantFunction
    value = 0
  []
[]

[LinearFVKernels]
  [mixture_ins_phase_2_time]
    type = LinearFVTimeDerivative
    variable = phase_2
  []
  [mixture_ins_phase_2_advection]
    type = LinearFVScalarAdvection
    advected_interp_method = upwind
    rhie_chow_user_object = ins_rhie_chow_interpolator
    u_slip = vel_slip_x
    v_slip = vel_slip_y
    variable = phase_2
  []
  [mixture_ins_phase_2_diffusion]
    type = LinearFVDiffusion
    diffusion_coeff = 1e-3
    use_nonorthogonal_correction = true
    variable = phase_2
  []
[]

[Materials]
  [mixture_phase_1_fraction]
    type = ParsedFunctorMaterial
    execute_on = ALWAYS
    expression = '1 - phase_2'
    functor_names = phase_2
    output_properties = phase_1
    outputs = none
    property_name = phase_1
  []
  [mixture_mixture_material]
    type = NSFVMixtureFunctorMaterial
    execute_on = ALWAYS
    limit_phase_fraction = true
    outputs = none
    phase_1_fraction = phase_2
    phase_1_names = '1.0 0.3 1 1'
    phase_2_names = '1e3 1.0 1 1'
    prop_names = 'rho_mixture mu_mixture cp_mixture k_mixture'
  []
  [mixture_slip_x]
    type = WCNSFV2PSlipVelocityFunctorMaterial
    execute_on = ALWAYS
    gravity = '0 -9.81 0'
    linear_coef_name = Darcy_coefficient
    momentum_component = x
    mu = mu_mixture
    outputs = none
    particle_diameter = 0.01
    rho = 1e3
    rho_d = 1.0
    slip_velocity_name = vel_slip_x
    u = vel_x
    v = vel_y
  []
  [mixture_slip_y]
    type = WCNSFV2PSlipVelocityFunctorMaterial
    execute_on = ALWAYS
    gravity = '0 -9.81 0'
    linear_coef_name = Darcy_coefficient
    momentum_component = y
    mu = mu_mixture
    outputs = none
    particle_diameter = 0.01
    rho = 1e3
    rho_d = 1.0
    slip_velocity_name = vel_slip_y
    u = vel_x
    v = vel_y
  []
  [mixture_dispersed_drag]
    type = NSFVDispersePhaseDragFunctorMaterial
    drag_coef_name = Darcy_coefficient
    execute_on = ALWAYS
    mu = mu_mixture
    outputs = none
    particle_diameter = 0.01
    rho = rho_mixture
    u = vel_x
    v = vel_y
  []
[]

[Variables]
  [phase_2]
    type = MooseLinearVariableFVReal
    family = MONOMIAL
    fv = true
    order = CONSTANT
    solver_sys = phi_system
  []
[]

[LinearFVBCs]
  [botttom-phase-2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'bottom'
    variable = phase_2
    functor = '0'
  []
  [top-phase-2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = phase_2
    functor = '1'
  []
[]

[Executioner]
  type = PIMPLE
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'

  end_time = 1e8
  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 2
    cutback_factor = 0.5
    dt = 1e-3
  []

  # Systems
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  active_scalar_systems = 'phi_system'
  momentum_equation_relaxation = 0.8
  active_scalar_equation_relaxation = '0.7'
  pressure_variable_relaxation = 0.3

  # We need to converge the problem to show conservation
  num_iterations = 200
  pressure_absolute_tolerance = 1e-10
  momentum_absolute_tolerance = 1e-10
  active_scalar_absolute_tolerance = '1e-10'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  active_scalar_petsc_options_iname = '-pc_type -pc_factor_shift_type' # -pc_hypre_type'
  active_scalar_petsc_options_value = 'lu NONZERO'
  momentum_l_abs_tol = 1e-13
  pressure_l_abs_tol = 1e-13
  active_scalar_l_abs_tol = 1e-13
  momentum_l_tol = 0
  pressure_l_tol = 0
  active_scalar_l_tol = 0
  # print_fields = true
  continue_on_max_its = true

  pin_pressure = true
  pressure_pin_value = 0.0
  pressure_pin_point = '0.0 0.0 0.0'
[]

[Outputs]
  exodus = false
  [out]
    type = CSV
    execute_on = 'FINAL'
  []
[]

[Postprocessors]
  [average_void]
    type = ElementAverageValue
    variable = 'phase_2'
  []
  [max_y_velocity]
    type = ElementExtremeValue
    variable = 'vel_y'
    value_type = max
  []
  [min_y_velocity]
    type = ElementExtremeValue
    variable = 'vel_y'
    value_type = min
  []
  [max_x_velocity]
    type = ElementExtremeValue
    variable = 'vel_x'
    value_type = max
  []
  [min_x_velocity]
    type = ElementExtremeValue
    variable = 'vel_x'
    value_type = min
  []
  [max_x_slip_velocity]
    type = ElementExtremeFunctorValue
    functor = 'vel_slip_x'
    value_type = max
  []
  [max_y_slip_velocity]
    type = ElementExtremeFunctorValue
    functor = 'vel_slip_y'
    value_type = max
  []
  [max_drag_coefficient_x]
    type = ElementExtremeFunctorValue
    functor = 'Darcy_coefficient_vec_out_x'
    value_type = max
  []
  [max_drag_coefficient_y]
    type = ElementExtremeFunctorValue
    functor = 'Darcy_coefficient_vec_out_y'
    value_type = max
  []
[]
