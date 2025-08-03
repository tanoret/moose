mu = 0.004
rho = 3260.0
advected_interp_method = 'upwind'

walls_fixed = 'left bottom right front back'
walls_all = '${walls_fixed} top'

### k-epsilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09
wall_treatment = 'two_layer'  # Options: eq_newton, eq_incremental, eq_linearized, neq

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    xmin = 0
    xmax = 1.0
    ymin = 0
    ymax = 1.0
    zmin = 0
    zmax = 1.0
    nx = 10
    ny = 10
    nz = 10
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system w_system pressure_system alpha_1_system TKE_system TKED_system'
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
    Re_y_star = 0.0 #60.0
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
    Re_y_star = 0.0 #60.0
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
    boundary = 'top'
    variable = vel_x
    functor = 1.0
  []
  [walls-v-impeller]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = vel_y
    functor = 0.0
  []
  [walls-w-impeller]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top'
    variable = vel_z
    functor = 0.0
  []

  [pressure-extrapolation]
    type = LinearFVExtrapolatedPressureBC
    boundary = '${walls_all}'
    variable = pressure
    use_two_term_expansion = true
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
    tke = 0.1 #TKE
    epsilon = 0.1 #TKED
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
