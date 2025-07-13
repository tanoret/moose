rho1 = 1000.0
rho2 = 1.0
advected_interp_method = 'upwind'
mu_1 = 1e-3
mu_2 = 1.48e-5
gravity = 9.81 #9.81

to_m = 0.146
domain_dims = ${fparse 4.0*to_m}
dam_dims_x = ${fparse 1*0.1461}
dam_dims_y = ${fparse 2*0.1461}

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${domain_dims}'
    dy = '${domain_dims}'
    ix = '500'
    iy = '500'
  []
[]

[Problem]
  linear_sys_names = 'u_system_p1 v_system_p1 u_system_p2 v_system_p2 pressure_system alpha_1_system alpha_2_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [rc_p1]
    type = RhieChowMassFluxMultiPhase
    u = vel_x_p1
    v = vel_y_p1
    pressure = pressure
    rho = ${rho1}
    p_diffusion_kernel = p_diffusion_p1
    check_executioner = false
    alpha = 'alpha_1'
    property_suffix = 'p1'
  []
  [rc_p2]
    type = RhieChowMassFluxMultiPhase
    u = vel_x_p2
    v = vel_y_p2
    pressure = pressure
    rho = ${rho2}
    p_diffusion_kernel = p_diffusion_p2
    check_executioner = false
    alpha = 'alpha_2'
    property_suffix = 'p2'
  []
[]

[Variables]
  [vel_x_p1]
    type = MooseLinearVariableFVReal
    solver_sys = u_system_p1
    initial_condition = 0.0
  []
  [vel_y_p1]
    type = MooseLinearVariableFVReal
    solver_sys = v_system_p1
    initial_condition = 0.0
  []
  [vel_x_p2]
    type = MooseLinearVariableFVReal
    solver_sys = u_system_p2
    initial_condition = 0.0
  []
  [vel_y_p2]
    type = MooseLinearVariableFVReal
    solver_sys = v_system_p2
    initial_condition = 0.0
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0.0
  []
  [alpha_1]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_1_system
  []
  [alpha_2]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_2_system
  []
[]

[LinearFVKernels]
  [u_time_p1]
    type = LinearFVTimeDerivative
    variable = vel_x_p1
    factor = ${rho1}
  []
  [v_time_p1]
    type = LinearFVTimeDerivative
    variable = vel_y_p1
    factor = ${rho1}
  []
  [u_advection_stress_p1]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_x_p1
    advected_interp_method = ${advected_interp_method}
    mu = ${mu_1}
    u = vel_x_p1
    v = vel_y_p1
    momentum_component = 'x'
    rhie_chow_user_object = 'rc_p1'
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [v_advection_stress_p1]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_y_p1
    advected_interp_method = ${advected_interp_method}
    mu = ${mu_1}
    u = vel_x_p1
    v = vel_y_p1
    momentum_component = 'y'
    rhie_chow_user_object = 'rc_p1'
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [u_pressure_p1]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_x_p1
    pressure = pressure
    alpha = 'alpha_1'
    momentum_component = 'x'
  []
  [v_pressure_p1]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_y_p1
    pressure = pressure
    alpha = 'alpha_1'
    momentum_component = 'y'
  []
  [v_gravity_p1]
    type = LinearFVSource
    variable = vel_y_p1
    source_density = 'alpha_1'
    scaling_factor = ${fparse -rho1 * gravity}
  []

  [u_time_p2]
    type = LinearFVTimeDerivative
    variable = vel_x_p2
    factor = ${rho2}
  []
  [v_time_p2]
    type = LinearFVTimeDerivative
    variable = vel_y_p2
    factor = ${rho2}
  []
  [u_advection_stress_p2]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_x_p2
    advected_interp_method = ${advected_interp_method}
    mu = ${mu_2}
    u = vel_x_p2
    v = vel_y_p2
    momentum_component = 'x'
    rhie_chow_user_object = 'rc_p2'
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [v_advection_stress_p2]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_y_p2
    advected_interp_method = ${advected_interp_method}
    mu = ${mu_2}
    u = vel_x_p2
    v = vel_y_p2
    momentum_component = 'y'
    rhie_chow_user_object = 'rc_p2'
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [u_pressure_p2]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_x_p2
    pressure = pressure
    alpha = 'alpha_2'
    momentum_component = 'x'
  []
  [v_pressure_p2]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_y_p2
    pressure = pressure
    alpha = 'alpha_2'
    momentum_component = 'y'
  []
  [v_gravity_p2]
    type = LinearFVSource
    variable = vel_y_p2
    source_density = 'alpha_2'
    scaling_factor = ${fparse -rho2 * gravity}
  []

  [p_diffusion_p1]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv_p1
    use_nonorthogonal_correction = false
  []
  [HbyA_divergence_p1]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA_p1
    force_boundary_execution = true
  []
  [p_diffusion_p2]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv_p2
    use_nonorthogonal_correction = false
  []
  [HbyA_divergence_p2]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA_p2
    force_boundary_execution = true
  []

  [alpha_1_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = alpha_1
    rho = 1.0 #${rho1}
    alpha = alpha_1
  []
  [alpha_1_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_1
    rhie_chow_user_object = 'rc_p1'
    c_alpha = 0.0
    rho = 1.0 #${rho1}
    u = vel_x_p1
    v = vel_x_p1
    u_mixture = vel_x_mixture
    v_mixture = vel_y_mixture
    limiter_method = 'vanLeer'
  []

  [alpha_2_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = alpha_2
    rho = 1.0 #${rho2}
    alpha = alpha_2
  []
  [alpha_2_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_2
    rhie_chow_user_object = 'rc_p2'
    c_alpha = 0.0
    rho = 1.0 #${rho2}
    u = vel_x_p2
    v = vel_x_p2
    u_mixture = vel_x_mixture
    v_mixture = vel_y_mixture
    limiter_method = 'vanLeer'
  []
[]

[LinearFVBCs]
  [walls-u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom'
    variable = vel_x_p1
    functor = 0.0
  []
  [walls-v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom'
    variable = vel_y_p1
    functor = 0.0
  []
  [outlet_u_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x_p1
    use_two_term_expansion = false
    boundary = 'top'
  []
  [outlet_v_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y_p1
    use_two_term_expansion = false
    boundary = 'top'
  []

  [walls-u_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom'
    variable = vel_x_p2
    functor = 0.0
  []
  [walls-v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left right bottom'
    variable = vel_y_p2
    functor = 0.0
  []
  [outlet_u_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x_p2
    use_two_term_expansion = false
    boundary = 'top'
  []
  [outlet_v_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y_p2
    use_two_term_expansion = false
    boundary = 'top'
  []

  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = pressure
    functor = 0.0
  []
  [pressure-extrapolation]
    type = LinearFVExtrapolatedPressureBC
    boundary = 'left right bottom'
    variable = pressure
    use_two_term_expansion = true
  []

  [walls_alpha_1]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_1
    functor = 0.0
    boundary = 'left right bottom'
  []
  [outlet_alpha_1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha_1
    use_two_term_expansion = false
    boundary = 'top'
  []

  [walls_alpha_2]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_2
    functor = 0.0
    boundary = 'left right bottom'
  []
  [outlet_alpha_2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha_2
    use_two_term_expansion = false
    boundary = 'top'
  []
[]

[ICs]
  [alpha_1]
    type = FunctionIC
    variable = 'alpha_1'
    function = alpha_1_init
  []
  [alpha_2]
    type = FunctionIC
    variable = 'alpha_2'
    function = alpha_2_init
  []
[]

[Functions]
  [alpha_1_init]
    type = ParsedFunction
    expression = 'if((x < ${dam_dims_x} & y < ${dam_dims_y}), 1.0, 0.0)'
  []
  [alpha_2_init]
    type = ParsedFunction
    expression = 'if((x > ${dam_dims_x} | y > ${dam_dims_y}), 1.0, 0.0)'
  []
[]

# [AuxVariables]
#   [alpha_2]
#     type = MooseLinearVariableFVReal
#   []
# []

# [AuxKernels]
#   [populate_alpha_2]
#     type = ParsedAux
#     variable = alpha_2
#     coupled_variables = 'alpha_1'
#     expression = 'min(max(1.0 - alpha_1, 0), 1)'
#     execute_on = 'NONLINEAR'
#   []
# []

[FunctorMaterials]
  [mixture_velocities]
    type = NSFVMixtureFunctorMaterial
    phase_1_names = 'vel_x_p1 vel_y_p2'
    phase_2_names = 'vel_x_p2 vel_y_p2'
    prop_names = 'vel_x_mixture vel_y_mixture'
    phase_1_fraction = 'alpha_1'
  []
[]

[Executioner]
  type = PIMPLEMultiPhase

  number_of_phases = 2

  momentum_l_abs_tol = 1e-12
  pressure_l_abs_tol = 1e-12
  momentum_l_tol = 1e-12
  pressure_l_tol = 1e-12

  rhie_chow_user_objects = 'rc_p1 rc_p2'
  momentum_systems = 'u_system_p1 v_system_p1; u_system_p2 v_system_p2'
  pressure_system = 'pressure_system'
  phase_systems = 'alpha_1_system alpha_2_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  phase_equation_relaxation = 0.5

  num_iterations = 100

  pressure_absolute_tolerance = 1e-11
  momentum_absolute_tolerance = 1e-11
  phase_absolute_tolerance = 1e-11

  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  phase_petsc_options_iname = '-pc_type -pc_hypre_type'
  phase_petsc_options_value = 'hypre boomeramg'

  print_fields = false
  continue_on_max_its = true
  dt = 0.0002
  num_steps = 500
  num_piso_iterations = 0
[]

[Outputs]
  exodus = true
[]
