rho1 = 1000.0
rho2 = 1.0

c_alpha = 0.0
advected_interp_method = 'upwind'
limiter_method = 'upwind' #'quick'

x_centroid = 50.0
y_centroid = 75.0
radius = 15.0

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '100.0'
    dy = '100.0'
    ix = '200'
    iy = '200'
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

[Functions]
  [u_vel_expression]
    type = ParsedFunction
    expression = '-sin(pi*x/100.0)^2 * sin(pi*y/50.0)'
  []
  [v_vel_expression]
    type = ParsedFunction
    expression = 'sin(pi*y/100.0)^2 * sin(pi*x/50.0)'
  []
[]

[LinearFVKernels]
  [u_reaction_p1]
    type = LinearFVReaction
    variable = vel_x_p1
  []
  [u_source_p1]
    type = LinearFVSource
    variable = vel_x_p1
    source_density = u_vel_expression
  []
  [v_reaction_p1]
    type = LinearFVReaction
    variable = vel_y_p1
  []
  [v_source_p1]
    type = LinearFVSource
    variable = vel_y_p1
    source_density = v_vel_expression
  []
  [u_reaction_p2]
    type = LinearFVReaction
    variable = vel_x_p2
  []
  [u_source_p2]
    type = LinearFVSource
    variable = vel_x_p2
    source_density = u_vel_expression
  []
  [v_reaction_p2]
    type = LinearFVReaction
    variable = vel_y_p2
  []
  [v_source_p2]
    type = LinearFVSource
    variable = vel_y_p2
    source_density = v_vel_expression
  []

  [p_diffusion_p1]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = one_mat
    use_nonorthogonal_correction = false
  []
  [p_diffusion_p2]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = one_mat
    use_nonorthogonal_correction = false
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
    c_alpha = ${c_alpha}
    rho = ${rho1}
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    activate_mules = false
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
    c_alpha = ${c_alpha}
    rho = ${rho2}
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    activate_mules = false
  []
[]

[LinearFVBCs]

  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = pressure
    functor = 0.0
  []

  # [walls_alpha_1]
  #   type = LinearFVAdvectionDiffusionFunctorNeumannBC
  #   variable = alpha_1
  #   functor = 0.0
  #   boundary = 'left right bottom'
  # []
  # [outlet_alpha_1]
  #   type = LinearFVAdvectionDiffusionOutflowBC
  #   variable = alpha_1
  #   use_two_term_expansion = false
  #   boundary = 'top'
  # []

  # [walls_alpha_2]
  #   type = LinearFVAdvectionDiffusionFunctorNeumannBC
  #   variable = alpha_2
  #   functor = 0.0
  #   boundary = 'left right bottom'
  # []
  # [outlet_alpha_2]
  #   type = LinearFVAdvectionDiffusionOutflowBC
  #   variable = alpha_2
  #   use_two_term_expansion = false
  #   boundary = 'top'
  # []
[]

[FunctorMaterials]
  [one_mat]
    type = GenericVectorFunctorMaterial
    prop_names = 'one_mat'
    prop_values = '1.0 1.0 1.0'
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
    expression = 'if(sqrt((x - ${x_centroid})^2 + (y - ${y_centroid})^2) < ${radius}, 0.0, 1.0)'
  []
  [alpha_2_init]
    type = ParsedFunction
    expression = 'if(sqrt((x - ${x_centroid})^2 + (y - ${y_centroid})^2) < ${radius}, 1.0, 0.0)'
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

  momentum_equation_relaxation = 1.0
  pressure_variable_relaxation = 0.3
  phase_equation_relaxation = 0.9

  num_iterations = 2

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
  dt = 0.001
  num_steps = 1000
  num_piso_iterations = 0

  # Interface tratment
  enforce_phase_sum = false
  activate_interface_shapening = false
  shapening_type = 'heaviside'
  smoothing_constant = 100.0
[]

[Outputs]
  exodus = true
[]
