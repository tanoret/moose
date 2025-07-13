rho1 = 100.0
rho2 = 10.0
advected_interp_method = 'upwind'
phase_1_in = 1.0
phase_2_in = 1.0
mu_1 = 0.2
mu_2 = 0.2

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1.'
    dy = '0.2'
    ix = '250'
    iy = '50'
  []
  [generate_inlet_top]
    type = ParsedGenerateSideset
    input = mesh
    fixed_normal = true
    normal = '-1 0 0'
    combinatorial_geometry = 'y >= 0.1 & x < 1e-5'
    new_sideset_name = 'inlet_top'
    include_only_external_sides = true
  []
  [generate_inlet_bottom]
    type = ParsedGenerateSideset
    input = generate_inlet_top
    fixed_normal = true
    normal = '-1 0 0'
    combinatorial_geometry = 'y <= 0.1 & x < 1e-5'
    new_sideset_name = 'inlet_bottom'
    include_only_external_sides = true
  []
  [delete_left]
    type = BoundaryDeletionGenerator
    input = generate_inlet_bottom
    boundary_names = left
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
    initial_condition = 0.5
  []
  [vel_y_p1]
    type = MooseLinearVariableFVReal
    solver_sys = v_system_p1
    initial_condition = 0.0
  []
  [vel_x_p2]
    type = MooseLinearVariableFVReal
    solver_sys = u_system_p2
    initial_condition = 0.5
  []
  [vel_y_p2]
    type = MooseLinearVariableFVReal
    solver_sys = v_system_p2
    initial_condition = 0.0
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0.2
  []
  [alpha_1]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_1_system
    initial_condition = 0.5
  []
  [alpha_2]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_2_system
    initial_condition = 0.5
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
    scaling_factor = ${fparse -rho1 * 9.81}
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
    scaling_factor = ${fparse -rho2 * 9.81}
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
    type = LinearFVTimeDerivative
    factor = ${rho1}
    variable = alpha_1
  []
  [alpha_1_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_1
    rhie_chow_user_object = 'rc_p1'
    c_alpha = 1e-4
    rho = ${rho1}
    u = vel_x_p1
    v = vel_x_p1
    u_mixture = vel_x_mixture
    v_mixture = vel_y_mixture
  []

  [alpha_2_time]
    type = LinearFVTimeDerivative
    factor = ${rho2}
    variable = alpha_2
  []
  [alpha_2_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_2
    rhie_chow_user_object = 'rc_p2'
    c_alpha = 1e-4
    rho = ${rho2}
    u = vel_x_p2
    v = vel_x_p2
    u_mixture = vel_x_mixture
    v_mixture = vel_y_mixture
  []
[]

[LinearFVBCs]
  [inlet-u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_top'
    variable = vel_x_p1
    functor = '1.1'
  []
  [inlet-v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_top'
    variable = vel_y_p1
    functor = '0.0'
  []
  [walls-u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom inlet_bottom'
    variable = vel_x_p1
    functor = 0.0
  []
  [walls-v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom inlet_bottom'
    variable = vel_y_p1
    functor = 0.0
  []
  [outlet_u_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x_p1
    use_two_term_expansion = false
    boundary = right
  []
  [outlet_v_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y_p1
    use_two_term_expansion = false
    boundary = right
  []

  [inlet-u_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_bottom'
    variable = vel_x_p2
    functor = '1.1'
  []
  [inlet-v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet_bottom'
    variable = vel_y_p2
    functor = '0.0'
  []
  [walls-u_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom inlet_top'
    variable = vel_x_p2
    functor = 0.0
  []
  [walls-v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom inlet_top'
    variable = vel_y_p2
    functor = 0.0
  []
  [outlet_u_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x_p2
    use_two_term_expansion = false
    boundary = right
  []
  [outlet_v_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y_p2
    use_two_term_expansion = false
    boundary = right
  []

  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'right'
    variable = pressure
    functor = 1.4
  []

  [inlet_alpha_1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = alpha_1
    functor = ${phase_1_in}
    boundary = 'inlet_top'
  []
  [walls_alpha_1]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_1
    functor = 0.0
    boundary = 'top bottom inlet_bottom'
  []
  [outlet_alpha_1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha_1
    use_two_term_expansion = false
    boundary = 'right'
  []

  [inlet_alpha_2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = alpha_2
    functor = ${phase_2_in}
    boundary = 'inlet_bottom'
  []
  [walls_alpha_2]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_2
    functor = 0.0
    boundary = 'top bottom inlet_top'
  []
  [outlet_alpha_2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha_2
    use_two_term_expansion = false
    boundary = 'right'
  []
[]

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

  num_iterations = 20

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
  dt = 1.0
  num_steps = 100
  num_piso_iterations = 0
[]

[Outputs]
  exodus = true
[]
