mu = 0.004
rho = 3260.0
advected_interp_method = 'upwind'

rpm = '1000'

walls_fixed = 'walls-bubbler walls-off-gas walls-spray walls'
walls_all = '${walls_fixed} walls-impeller'

### k-epsilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09
wall_treatment = 'two_layer'  # Options: eq_newton, eq_incremental, eq_linearized, neq

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = 'mesh_pump_bowl_test.e'
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system TKE_system TKED_system alpha_1_system'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc]
    type = RhieChowMassFluxMultiPhase
    u = vel_x
    v = vel_y
    w = vel_z
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
    check_executioner = false
    alpha = '1.0'
    # property_suffix = 'p1'
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    initial_condition = 0.0
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    solver_sys = v_system
    initial_condition = 0.0
  []
  [vel_z]
    type = MooseLinearVariableFVReal
    solver_sys = w_system
    initial_condition = 0.0
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0.2
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = 0.1
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = 0.1
  []
  [alpha_1]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_1_system
  []
[]

[LinearFVKernels]
  [u_time]
    type = LinearFVTimeDerivative
    variable = vel_x
    factor = ${rho}
  []
  [v_time]
    type = LinearFVTimeDerivative
    variable = vel_y
    factor = ${rho}
  []
  [w_time]
    type = LinearFVTimeDerivative
    variable = vel_z
    factor = ${rho}
  []
  [u_advection_stress]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [v_advection_stress]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [w_advection_stress]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_z
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    w = vel_z
    momentum_component = 'z'
    rhie_chow_user_object = 'rc'
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [w_diffusion]
    type = LinearFVDiffusion
    variable = vel_z
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [u_pressure]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_x
    pressure = pressure
    alpha = 1.0
    momentum_component = 'x'
  []
  [v_pressure]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_y
    pressure = pressure
    alpha = 1.0
    momentum_component = 'y'
  []
  [w_pressure]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_z
    pressure = pressure
    alpha = 1.0
    momentum_component = 'z'
  []

  [p_diffusion]
    type = LinearFVAnisotropicDiffusion
    variable = pressure
    diffusion_tensor = Ainv
    use_nonorthogonal_correction = false
  []
  [HbyA_divergence]
    type = LinearFVDivergence
    variable = pressure
    face_flux = HbyA
    force_boundary_execution = true
  []

  [TKE_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = TKE
    rho = ${rho}
    alpha = 1.0
  []
  [TKE_advection]
    type = LinearFVTurbulentMultiPhaseAdvection
    variable = TKE
    alpha = 1.0
  []
  [TKE_diffusion]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    alpha = 1.0
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink]
    type = LinearFVTKEMultiPhaseSourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    w = vel_z
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    wall_distance = 'd'
    alpha = 1.0
  []

  [TKED_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = TKED
    rho = ${rho}
    alpha = 1.0
  []
  [TKED_advection]
    type = LinearFVTurbulentMultiPhaseAdvection
    variable = TKED
    tke = 'TKE'
    rho = ${rho}
    mu = ${mu}
    alpha = 1.0
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_diffusion]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKED
    diffusion_coeff = '${mu}'
    use_nonorthogonal_correction = false
    tke = 'TKE'
    rho = ${rho}
    mu = 0.004 #'${mu}'
    alpha = 1.0
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    tke = 'TKE'
    rho = ${rho}
    mu = ${mu}
    alpha = 1.0
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_source_sink]
    type = LinearFVTKEDMultiPhaseSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    w = vel_z
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    wall_distance = 'd'
    alpha = 1.0
  []

  [alpha_1_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = alpha_1
    rho = ${rho}
    alpha = 1.0
  []
  # [alpha_1_advection]
  #   type = LinearFVMultiPhaseFractionAdvection
  #   variable = alpha_1
  #   rhie_chow_user_object = 'rc_p1'
  #   c_alpha = ${c_alpha}
  #   rho = ${rho}
  #   advected_interp_method = ${advected_interp_method}
  #   limiter_method = ${limiter_method}
  #   use_nonorthogonal_correction = false
  # []
[]

[LinearFVBCs]
  # [inlet-u]
  #   type = LinearFVAdvectionDiffusionFunctorDirichletBC
  #   boundary = 'inlet'
  #   variable = vel_x
  #   functor = '0.0'
  # []
  # [inlet-v]
  #   type = LinearFVAdvectionDiffusionFunctorDirichletBC
  #   boundary = 'inlet'
  #   variable = vel_y
  #   functor = '0.0'
  # []
  # [inlet-w]
  #   type = LinearFVAdvectionDiffusionFunctorDirichletBC
  #   boundary = 'inlet'
  #   variable = vel_z
  #   functor = '1.0'
  # []

  [walls-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_x
    functor = 0.0
  []
  [walls-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_y
    functor = 0.0
  []
  [walls-w]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_z
    functor = 0.0
  []

  [walls-u-impeller]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_x
    functor = 'vel_x_impeller'
  []
  [walls-v-impeller]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_y
    functor = 'vel_y_impeller'
  []
  [walls-w-impeller]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_z
    functor = 0.0
  []

  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet outlet'
    variable = pressure
    functor = 1.4
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x
    use_two_term_expansion = false
    boundary = 'inlet outlet'
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y
    use_two_term_expansion = false
    boundary = 'inlet outlet'
  []
  [outlet_w]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_z
    use_two_term_expansion = false
    boundary = 'inlet outlet'
  []

  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'inlet outlet'
    variable = TKE
    use_two_term_expansion = false
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'inlet outlet'
    variable = TKED
    use_two_term_expansion = false
  []
  [outlet_alpha_1]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'inlet outlet'
    variable = alpha_1
    use_two_term_expansion = false
  []
[]

[AuxVariables]
  [mu_t]
    type = MooseLinearVariableFVReal
    initial_condition = '0.01'
  []
  # [yplus]
  #   type = MooseLinearVariableFVReal
  # []
  [d]
    type = MooseLinearVariableFVReal
    initial_condition = '1.0'
  []
[]

[AuxKernels]
  [compute_mu_t]
    type = kEpsilonViscosityAux
    variable = mu_t
    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    w = vel_z
    walls = '${walls_all}'
    wall_treatment = ${wall_treatment}
    wall_distance = 'd'
    execute_on = 'NONLINEAR'
  []
  # [compute_y_plus]
  #   type = RANSYPlusAux
  #   variable = yplus
  #   tke = TKE
  #   mu = ${mu}
  #   rho = ${rho}
  #   u = vel_x
  #   v = vel_y
  #   w = vel_z
  #   walls = '${walls_all}'
  #   wall_treatment = 'eq_newton'
  #   execute_on = 'NONLINEAR'
  # []
  [compute_wall_distance]
    type = WallDistanceAux
    variable = d
    walls = '${walls_all}'
    execute_on = 'INITIAL'
  []
[]

[Functions]
  [radius]
    type = ParsedFunction
    expression = 'sqrt(x^2 + y^2)'
  []
  [vel_x_impeller]
    type = ParsedFunction
    expression = '-${rpm} * r * 2*pi/60 * y/r'
    symbol_names = 'r'
    symbol_values = 'radius'
  []
  [vel_y_impeller]
    type = ParsedFunction
    expression = '${rpm} * r * 2*pi/60 * x/r'
    symbol_names = 'r'
    symbol_values = 'radius'
  []
[]

[Executioner]
  type = PIMPLEMultiPhase

  number_of_phases = 1

  momentum_l_abs_tol = 1e-10
  pressure_l_abs_tol = 1e-10
  turbulence_l_abs_tol = 1e-14
  phase_l_abs_tol = 1e-14

  momentum_l_tol = 0
  pressure_l_tol = 0
  turbulence_l_tol = 0
  phase_l_tol = 0

  rhie_chow_user_objects = 'rc'
  momentum_systems = 'u_system v_system w_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system'
  phase_systems = 'alpha_1_system'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.25 0.25'
  turbulence_field_relaxation = '0.25 0.25'
  phase_equation_relaxation = 0.9

  num_iterations = 10
  dt = 0.1
  num_steps = 50
  num_piso_iterations = 0

  pressure_absolute_tolerance = 1e-10
  momentum_absolute_tolerance = 1e-10
  turbulence_absolute_tolerance = '1e-12 1e-12'
  phase_absolute_tolerance = 1e-11

  momentum_petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type'
  momentum_petsc_options_value = 'hypre boomeramg 4 1 0.1 0.6 HMIS ext+i'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type'
  pressure_petsc_options_value = 'hypre boomeramg 2 1 0.1 0.6 HMIS ext+i'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type'
  turbulence_petsc_options_value = 'hypre boomeramg 2 1 0.1 0.6 HMIS ext+i'
  phase_petsc_options_iname = '-pc_type -pc_hypre_type -pc_hypre_boomeramg_agg_nl -pc_hypre_boomeramg_agg_num_paths -pc_hypre_boomeramg_truncfactor -pc_hypre_boomeramg_strong_threshold -pc_hypre_boomeramg_coarsen_type -pc_hypre_boomeramg_interp_type'
  phase_petsc_options_value = 'hypre boomeramg 2 1 0.1 0.6 HMIS ext+i'

  print_fields = false
  continue_on_max_its = true

[]

[Outputs]
  exodus = true
[]
