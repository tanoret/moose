# phase 1 is water, phase 2 is gas
# Hosokawa and Tomiyama (2013)

# flow regime bubbly
# void fraction(%)       0.18     0.71     1.56     1.27
# liquid velocity(m/s)   25       6.25     2.88     3.54
# gas velocity(m/s)      2.22e-1  2.22e-1  2.31e-1  2.36e-1
# bubble diameter(m)     3.48e-3  3.52e-3  3.59e-3 2.62e-3

mu_1 = 8.899e-4
mu_2 = 1.831e-5
rho_1 = 998.19
rho_2 = 1.185
pipe_diameter = 0.02 # m
length = 2 # m
U_1 = 2.88   # m/s
U_2 = 2.31e-1   # m/s
g = -9.81
inlet_phase_2 = 0.0156
advected_interp_method = 'upwind'
velocity_interp_method = 'rc'
average_bubble_diameter = 3.59e-3

rho1 = 100.0
rho2 = 100.0
advected_interp_method = 'upwind'
phase_1_in = 0.7
phase_2_in = 0.3
mu_1 = 2.6
mu_2 = 2.6

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '1.'
    dy = '0.2'
    ix = '100'
    iy = '20'
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
    initial_condition = ${phase_1_in}
  []
  [alpha_2]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_2_system
    initial_condition = ${phase_2_in}
  []
[]

[LinearFVKernels]
  # [u_time_p1]
  #   type = LinearFVTimeDerivative
  #   variable = vel_x_p1
  #   factor = ${rho1}
  # []
  # [v_time_p1]
  #   type = LinearFVTimeDerivative
  #   variable = vel_y_p1
  #   factor = ${rho1}
  # []
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

  # [u_time_p2]
  #   type = LinearFVTimeDerivative
  #   variable = vel_x_p2
  #   factor = ${rho2}
  # []
  # [v_time_p2]
  #   type = LinearFVTimeDerivative
  #   variable = vel_y_p2
  #   factor = ${rho2}
  # []
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

  # [alpha_1_time]
  #   type = LinearFVTimeDerivative
  #   factor = ${rho1}
  #   variable = alpha_1
  # []
  [alpha_1_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_1
    rhie_chow_user_object = 'rc_p1'
  []

  # [alpha_2_time]
  #   type = LinearFVTimeDerivative
  #   factor = ${rho2}
  #   variable = alpha_2
  # []
  [alpha_2_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_2
    rhie_chow_user_object = 'rc_p2'
  []
[]

[LinearFVBCs]
  [inlet-u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_x_p1
    functor = '1.1'
  []
  [inlet-v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_y_p1
    functor = '0.0'
  []
  [walls-u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_x_p1
    functor = 0.0
  []
  [walls-v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
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
    boundary = 'left'
    variable = vel_x_p2
    functor = '1.1'
  []
  [inlet-v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_y_p2
    functor = '0.0'
  []
  [walls-u_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_x_p2
    functor = 0.0
  []
  [walls-v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
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
    boundary = 'left'
  []
  [walls_alpha_1]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_1
    functor = 0.0
    boundary = 'top bottom'
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
    boundary = 'left'
  []
  [walls_alpha_2]
    type = LinearFVAdvectionDiffusionFunctorNeumannBC
    variable = alpha_2
    functor = 0.0
    boundary = 'top bottom'
  []
  [outlet_alpha_2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = alpha_2
    use_two_term_expansion = false
    boundary = 'right'
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

  momentum_equation_relaxation = 0.8
  pressure_variable_relaxation = 0.3
  phase_equation_relaxation = 0.7

  num_iterations = 200

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
  num_steps = 1
  num_piso_iterations = 0
[]

[Postprocessors]
  [Re]
    type = ParsedPostprocessor
    expression = '${rho_1} * ${fparse pipe_diameter} * ${U_1}'
    pp_names = 'Re'
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
  exodus = true
[]
