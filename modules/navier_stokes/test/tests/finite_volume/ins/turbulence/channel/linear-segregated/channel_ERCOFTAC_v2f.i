H = 1 #halfwidth of the channel
L = 100

Re = 13700

rho = 1
bulk_u = 1
mu = '${fparse rho * bulk_u * 2 * H / Re}'

advected_interp_method = 'upwind'

### k-epsilon Closure Parameters ###
sigma_k = 1.0
sigma_theta_2 = 1.0
sigma_eps = 1.3
C1_eps = 1.4
C2_eps = 1.9
C_mu = 0.09
C_mu_theta_2 = 0.22
CL = 0.23
C_eta = 70.0

### Initial and Boundary Conditions ###
intensity = 0.01
k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / (2*H)}'

### Modeling parameters ###
# bulk_wall_treatment = false
walls = 'top bottom'
wall_treatment = 'eq_newton'  # Options: eq_newton, eq_incremental, eq_linearized, neq

[Mesh]
  [block_1]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = 0
    ymax = ${H}
    nx = 100
    ny = 10
    bias_y = 0.8
  []
  [block_2]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = ${L}
    ymin = ${fparse -H}
    ymax = 0
    nx = 100
    ny = 10
    bias_y = ${fparse 1/0.8}
  []
  [smg]
    type = StitchedMeshGenerator
    inputs = 'block_1 block_2'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'bottom top'
    merge_boundaries_with_same_name = true
  []
[]

[Problem]
  linear_sys_names = 'u_system v_system pressure_system TKE_system TKED_system v2_system f_system'
  previous_nl_solution_required = true
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc]
    type = RhieChowMassFlux
    u = vel_x
    v = vel_y
    pressure = pressure
    rho = ${rho}
    p_diffusion_kernel = p_diffusion
  []
[]

[Variables]
  [vel_x]
    type = MooseLinearVariableFVReal
    initial_condition = ${bulk_u}
    solver_sys = u_system
  []
  [vel_y]
    type = MooseLinearVariableFVReal
    initial_condition = 0
    solver_sys = v_system
  []
  [pressure]
    type = MooseLinearVariableFVReal
    initial_condition = 1e-8
    solver_sys = pressure_system
  []
  [TKE]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system
    initial_condition = ${k_init}
  []
  [TKED]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system
    initial_condition = ${eps_init}
  []
  [theta_squared]
    type = MooseLinearVariableFVReal
    solver_sys = v2_system
    initial_condition = ${k_init}
  []
  [f]
    type = MooseLinearVariableFVReal
    solver_sys = f_system
    initial_condition = 1.0
  []
[]

[LinearFVKernels]
  [u_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'x'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
  []
  [u_diffusion]
    type = LinearFVDiffusion
    variable = vel_x
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [u_pressure]
    type = LinearFVMomentumPressure
    variable = vel_x
    pressure = pressure
    momentum_component = 'x'
  []
  [v_advection_stress]
    type = LinearWCNSFVMomentumFlux
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t'
    u = vel_x
    v = vel_y
    momentum_component = 'y'
    rhie_chow_user_object = 'rc'
    use_nonorthogonal_correction = false
    use_deviatoric_terms = yes
  []
  [v_diffusion]
    type = LinearFVDiffusion
    variable = vel_y
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [v_pressure]
    type = LinearFVMomentumPressure
    variable = vel_y
    pressure = pressure
    momentum_component = 'y'
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

  [TKE_advection]
    type = LinearFVTurbulentAdvection
    variable = TKE
  []
  [TKE_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKE
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_k}
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink]
    type = LinearFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C_pl = 1e10
    walls = ${walls}
    wall_treatment = ${wall_treatment}
  []

  [TKED_advection]
    type = LinearFVTurbulentAdvection
    variable = TKED
    walls = ${walls}
  []
  [TKED_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = ${mu}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = TKED
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    walls = ${walls}
  []
  [TKED_source_sink]
    type = LinearFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    tke = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    C_pl = 1e10
    v2f_bool = true
    theta_squared = 'theta_squared'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
  []

  [theta_squared_advection]
    type = LinearFVTurbulentAdvection
    variable = theta_squared
    walls = ${walls}
  []
  [theta_squared_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = theta_squared
    diffusion_coeff = ${mu}
    walls = ${walls}
    use_nonorthogonal_correction = false
  []
  [theta_squared_turb_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = theta_squared
    diffusion_coeff = 'mu_t'
    scaling_coeff = ${sigma_theta_2}
    walls = ${walls}
    use_nonorthogonal_correction = false
  []
  [theta_squared_source_sink]
    type = LinearFVThetaSquaredSourceSink
    variable = theta_squared
    u = vel_x
    v = vel_y
    tke = TKE
    epsilon = TKED
    f = f
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
  []

  [f_diffusion]
    type = LinearFVTurbulentDiffusion
    variable = f
    diffusion_coeff = 'L2'
    use_nonorthogonal_correction = true
    walls = ${walls}
  []
  [f_source_sink]
    type = LinearFVEllipticBlendingSourceSink
    variable = f
    u = vel_x
    v = vel_y
    tke = TKE
    epsilon = TKED
    theta_squared = theta_squared
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
  []
[]

[LinearFVBCs]
  [inlet-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_x
    functor = '${bulk_u}'
  []
  [inlet-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = vel_y
    functor = '0.0'
  []
  [walls-u]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_x
    functor = 0.0
  []
  [walls-v]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'top bottom'
    variable = vel_y
    functor = 0.0
  []
  [outlet_u]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = vel_x
    use_two_term_expansion = false
  []
  [outlet_v]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = vel_y
    use_two_term_expansion = false
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'right'
    variable = pressure
    functor = 0.0
  []

  [inlet_TKE]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = TKE
    functor = '${k_init}'
  []
  [outlet_TKE]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = TKE
    use_two_term_expansion = false
  []
  [inlet_theta_squared]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = theta_squared
    functor = '${k_init}'
  []
  [outlet_theta_squared]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = theta_squared
    use_two_term_expansion = false
  []
  [inlet_TKED]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'left'
    variable = TKED
    functor = '${eps_init}'
  []
  [outlet_TKED]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'right'
    variable = TKED
    use_two_term_expansion = false
  []
[]

[FVBCs]
  [walls_mu_t]
    type = INSFVTurbulentViscosityWallFunction
    boundary = 'bottom top'
    variable = mu_t
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t'
    tke = TKE
    wall_treatment = ${wall_treatment}
  []
[]

[AuxVariables]
  [time_scale]
    type = MooseVariableFVReal
    initial_condition = '${fparse ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
  [time_scale_realizable]
    type = MooseVariableFVReal
    initial_condition = '${fparse ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
  [mu_t_keps]
    type = MooseVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
  [mu_t]
    type = MooseVariableFVReal
    initial_condition = '${fparse rho * C_mu * ${k_init}^2 / eps_init}'
    two_term_boundary_expansion = false
  []
  [mixing_length_viscosity]
    type = MooseVariableFVReal
    initial_condition = '1.0'
    two_term_boundary_expansion = false
  []
  [strain_rate]
    type = MooseVariableFVReal
    initial_condition = '1.0'
    two_term_boundary_expansion = false
  []
  [realizable_scale]
    type = MooseVariableFVReal
    initial_condition = '1.0'
    two_term_boundary_expansion = false
  []
  [yplus]
    type = MooseVariableFVReal
    two_term_boundary_expansion = false
  []
  [L2]
    type = MooseVariableFVReal
    initial_condition = '0.1'
    two_term_boundary_expansion = false
  []
  [u_var]
    type = INSFVVelocityVariable
    initial_condition = '1.0'
    two_term_boundary_expansion = false
  []
  [v_var]
    type = INSFVVelocityVariable
    initial_condition = '1.0'
    two_term_boundary_expansion = false
  []
[]

[AuxKernels]
  [compute_u_var]
    type = FunctorAux
    functor = 'vel_x'
    variable = 'u_var'
    execute_on = 'NONLINEAR'
  []
  [compute_v_var]
    type = FunctorAux
    functor = 'vel_y'
    variable = 'v_var'
    execute_on = 'NONLINEAR'
  []
  [compute_mixing_length_viscosity]
    type = INSFVMixingLengthTurbulentViscosityAux
    variable = mixing_length_viscosity
    mixing_length = 1.0
    u = u_var
    v = v_var
    execute_on = 'NONLINEAR'
  []
  [compute_strain_rate]
    type = ParsedAux
    variable = strain_rate
    functor_names = 'mixing_length_viscosity'
    functor_symbols = 'mixing_length_viscosity'
    expression = 'sqrt(mixing_length_viscosity)'
    execute_on = 'NONLINEAR'
  []
  [compute_realizable_scale]
    type = ParsedAux
    variable = realizable_scale
    functor_names = 'TKE theta_squared strain_rate'
    functor_symbols = 'TKE theta_squared strain_rate'
    expression = '0.6*TKE/(sqrt(3)*${C_mu_theta_2}*theta_squared*strain_rate)'
    execute_on = 'NONLINEAR'
  []
  [compute_time_scale_realizable]
    type = ParsedAux
    variable = time_scale_realizable
    functor_names = 'TKE TKED realizable_scale'
    functor_symbols = 'TKE TKED realizable_scale'
    expression = 'max(min(TKE/TKED, realizable_scale), 6.0*sqrt(${mu}/${rho}/TKED))'
    execute_on = 'NONLINEAR'
  []
  [compute_time_scale]
    type = ParsedAux
    variable = time_scale
    functor_names = 'TKE TKED'
    functor_symbols = 'TKE TKED'
    expression = 'max(TKE/TKED, 6.0*sqrt(${mu}/${rho}/TKED))'
    execute_on = 'NONLINEAR'
  []
  [compute_mu_t]
    type = ParsedAux
    variable = mu_t
    functor_names = 'TKE TKED theta_squared time_scale time_scale_realizable mu_t_keps'
    functor_symbols = 'TKE TKED theta_squared time_scale time_scale_realizable mu_t_keps'
    expression = '${rho} * min(${C_mu}*TKE*time_scale_realizable, ${C_mu_theta_2}*theta_squared*time_scale)'
    # expression = 'min(mu_t_keps, ${rho}*${C_mu_theta_2}*theta_squared*time_scale)'
    execute_on = 'NONLINEAR'
  []
  [compute_mu_t_keps]
    type = kEpsilonViscosityAux
    variable = mu_t_keps
    C_mu = ${C_mu}
    tke = TKE
    epsilon = TKED
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    bulk_wall_treatment = false
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
    mu_t_ratio_max = 1e20
  []
  [compute_y_plus]
    type = RANSYPlusAux
    variable = yplus
    tke = TKE
    mu = ${mu}
    rho = ${rho}
    u = vel_x
    v = vel_y
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    execute_on = 'NONLINEAR'
  []
  [compute_L2]
    type = ParsedAux
    variable = L2
    functor_names = 'TKE TKED time_scale_realizable'
    functor_symbols = 'TKE TKED time_scale_realizable'
    expression = '(${CL} * max(min(TKE^1.5/TKED, time_scale_realizable*TKE^0.5/0.6), ${C_eta}*((${mu}/${rho})^3/TKED)^0.25))^2'
    # expression = '(${CL}*max(TKE^1.5/TKED,${C_eta}*((${mu}/${rho})^3/TKED)^0.25))^2'
    execute_on = 'NONLINEAR'
  []
[]

[Executioner]
  type = SIMPLE

  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system TKED_system f_system v2_system'

  momentum_l_abs_tol = 1e-14
  pressure_l_abs_tol = 1e-14
  turbulence_l_abs_tol = 1e-14
  momentum_l_tol = 1e-14
  pressure_l_tol = 1e-14
  turbulence_l_tol = 1e-14

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.25 0.25 0.25 0.25'
  num_iterations = 1000
  pressure_absolute_tolerance = 1e-12
  momentum_absolute_tolerance = 1e-12
  turbulence_absolute_tolerance = '1e-12 1e-12 1e-12 1e-12'
  momentum_petsc_options_iname = '-pc_type -pc_hypre_type'
  momentum_petsc_options_value = 'hypre boomeramg'
  pressure_petsc_options_iname = '-pc_type -pc_hypre_type'
  pressure_petsc_options_value = 'hypre boomeramg'
  turbulence_petsc_options_iname = '-pc_type -pc_hypre_type'
  turbulence_petsc_options_value = 'hypre boomeramg'

  print_fields = false
  continue_on_max_its = true
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
