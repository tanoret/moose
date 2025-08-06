rho1 = 3260.0
rho2 = 1.0
mu1 = 0.004
mu2 = 1.48e-5
gravity = 9.81

c_alpha = 0.1
advected_interp_method = 'upwind'
limiter_method = 'vanLeer'

MULES_iterations = 1

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

phase_2_inlet = 0.01
phase_1_inlet = ${fparse 1.0 - phase_2_inlet}

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = 'mesh_pump_bowl_test.e'
  []
[]

[Problem]
  linear_sys_names = 'u_system_p1 v_system_p1 w_system_p1
                      u_system_p2 v_system_p2 w_system_p2
                      pressure_system
                      TKE_system_p1 TKED_system_p1
                      TKE_system_p2 TKED_system_p2
                      alpha_system_p1 alpha_system_p2'
  previous_nl_solution_required = true
[]

[GlobalParams]
  advected_interp_method = ${advected_interp_method}
[]

[UserObjects]
  [rc_p1]
    type = RhieChowMassFluxMultiPhase
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
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
    w = vel_z_p2
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
  [vel_z_p1]
    type = MooseLinearVariableFVReal
    solver_sys = w_system_p1
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
  [vel_z_p2]
    type = MooseLinearVariableFVReal
    solver_sys = w_system_p2
    initial_condition = 0.0
  []
  [pressure]
    type = MooseLinearVariableFVReal
    solver_sys = pressure_system
    initial_condition = 0.2
  []
  [TKE_p1]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system_p1
    initial_condition = 0.1
  []
  [TKED_p1]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system_p1
    initial_condition = 0.1
  []
  [TKE_p2]
    type = MooseLinearVariableFVReal
    solver_sys = TKE_system_p2
    initial_condition = 0.1
  []
  [TKED_p2]
    type = MooseLinearVariableFVReal
    solver_sys = TKED_system_p2
    initial_condition = 0.1
  []
  [alpha_1]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_system_p1
    initial_condition = '${phase_1_inlet}'
  []
  [alpha_2]
    type = MooseLinearVariableFVReal
    solver_sys = alpha_system_p2
    initial_condition = '${phase_2_inlet}'
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
  [w_time_p1]
    type = LinearFVTimeDerivative
    variable = vel_z_p1
    factor = ${rho1}
  []
  [u_advection_stress_p1]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_x_p1
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t_p1'
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
    momentum_component = 'x'
    rhie_chow_user_object = 'rc_p1'
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [v_advection_stress_p1]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_y_p1
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t_p1'
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
    momentum_component = 'y'
    rhie_chow_user_object = 'rc_p1'
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [w_advection_stress_p1]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_z_p1
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t_p1'
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
    momentum_component = 'z'
    rhie_chow_user_object = 'rc_p1'
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [u_diffusion_p1]
    type = LinearFVDiffusion
    variable = vel_x_p1
    diffusion_coeff = ${mu1}
    use_nonorthogonal_correction = false
  []
  [v_diffusion_p1]
    type = LinearFVDiffusion
    variable = vel_y_p1
    diffusion_coeff = ${mu1}
    use_nonorthogonal_correction = false
  []
  [w_diffusion_p1]
    type = LinearFVDiffusion
    variable = vel_z_p1
    diffusion_coeff = ${mu1}
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
  [w_pressure_p1]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_z_p1
    pressure = pressure
    alpha = 'alpha_1'
    momentum_component = 'z'
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
  [w_time_p2]
    type = LinearFVTimeDerivative
    variable = vel_z_p2
    factor = ${rho2}
  []
  [u_advection_stress_p2]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_x_p2
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t_p2'
    u = vel_x_p2
    v = vel_y_p2
    w = vel_z_p2
    momentum_component = 'x'
    rhie_chow_user_object = 'rc_p2'
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [v_advection_stress_p2]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_y_p2
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t_p2'
    u = vel_x_p2
    v = vel_y_p2
    w = vel_z_p2
    momentum_component = 'y'
    rhie_chow_user_object = 'rc_p2'
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [w_advection_stress_p2]
    type = LinearWCNSFVMultiPhaseMomentumFlux
    variable = vel_z_p2
    advected_interp_method = ${advected_interp_method}
    mu = 'mu_t_p2'
    u = vel_x_p2
    v = vel_y_p2
    w = vel_z_p2
    momentum_component = 'z'
    rhie_chow_user_object = 'rc_p2'
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [u_diffusion_p2]
    type = LinearFVDiffusion
    variable = vel_x_p2
    diffusion_coeff = ${mu2}
    use_nonorthogonal_correction = false
  []
  [v_diffusion_p2]
    type = LinearFVDiffusion
    variable = vel_y_p2
    diffusion_coeff = ${mu2}
    use_nonorthogonal_correction = false
  []
  [w_diffusion_p2]
    type = LinearFVDiffusion
    variable = vel_z_p2
    diffusion_coeff = ${mu2}
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
  [w_pressure_p2]
    type = LinearFVMultiPhaseMomentumPressure
    variable = vel_z_p2
    pressure = pressure
    alpha = 'alpha_2'
    momentum_component = 'z'
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

  [TKE_time_p1]
    type = LinearFVMultiPhaseTimeDerivative
    variable = TKE_p1
    rho = ${rho1}
    alpha = 'alpha_1'
  []
  [TKE_advection_p1]
    type = LinearFVTurbulentMultiPhaseAdvection
    variable = TKE_p1
    rhie_chow_user_object = 'rc_p1'
    alpha = 'alpha_1'
  []
  [TKE_diffusion_p1]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKE_p1
    diffusion_coeff = ${mu1}
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion_p1]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKE_p1
    diffusion_coeff = 'mu_t_p1'
    scaling_coeff = ${sigma_k}
    alpha = 'alpha_1'
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink_p1]
    type = LinearFVTKEMultiPhaseSourceSink
    variable = TKE_p1
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
    epsilon = TKED_p1
    rho = ${rho1}
    mu = ${mu1}
    wall_distance = 'd'
    alpha = 'alpha_1'
  []

  [TKED_time_p1]
    type = LinearFVMultiPhaseTimeDerivative
    variable = TKED_p1
    rho = ${rho1}
    alpha = 'alpha_1'
  []
  [TKED_advection_p1]
    type = LinearFVTurbulentMultiPhaseAdvection
    variable = TKED_p1
    rhie_chow_user_object = 'rc_p1'
    tke = 'TKE_p1_aux'
    rho = ${rho1}
    mu = ${mu1}
    alpha = 'alpha_1'
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_diffusion_p1]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKED_p1
    diffusion_coeff = '${mu1}'
    use_nonorthogonal_correction = false
    tke = 'TKE_p1_aux'
    rho = ${rho1}
    mu = 0.004 #'${mu}'
    alpha = 'alpha_1'
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_turb_diffusion_p1]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKED_p1
    diffusion_coeff = 'mu_t_p1'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    tke = 'TKE_p1_aux'
    rho = ${rho1}
    mu = ${mu1}
    alpha = 'alpha_1'
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_source_sink_p1]
    type = LinearFVTKEDMultiPhaseSourceSink
    variable = TKED_p1
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
    tke = TKE_p1
    rho = ${rho1}
    mu = ${mu1}
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    wall_distance = 'd'
    alpha = 'alpha_1'
  []

  [TKE_time_p2]
    type = LinearFVMultiPhaseTimeDerivative
    variable = TKE_p2
    rho = ${rho2}
    alpha = 'alpha_2'
  []
  [TKE_advection_p2]
    type = LinearFVTurbulentMultiPhaseAdvection
    variable = TKE_p2
    rhie_chow_user_object = 'rc_p2'
    alpha = 'alpha_2'
  []
  [TKE_diffusion_p2]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKE_p2
    diffusion_coeff = ${mu2}
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [TKE_turb_diffusion_p2]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKE_p2
    diffusion_coeff = 'mu_t_p2'
    scaling_coeff = ${sigma_k}
    alpha = 'alpha_2'
    use_nonorthogonal_correction = false
  []
  [TKE_source_sink_p2]
    type = LinearFVTKEMultiPhaseSourceSink
    variable = TKE_p2
    u = vel_x_p2
    v = vel_y_p2
    w = vel_z_p2
    epsilon = TKED_p2
    rho = ${rho2}
    mu = ${mu2}
    wall_distance = 'd'
    alpha = 'alpha_2'
  []

  [TKED_time_p2]
    type = LinearFVMultiPhaseTimeDerivative
    variable = TKED_p2
    rho = ${rho2}
    alpha = 'alpha_2'
  []
  [TKED_advection_p2]
    type = LinearFVTurbulentMultiPhaseAdvection
    variable = TKED_p2
    rhie_chow_user_object = 'rc_p2'
    tke = 'TKE_p2_aux'
    rho = ${rho2}
    mu = ${mu2}
    alpha = 'alpha_2'
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_diffusion_p2]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKED_p2
    diffusion_coeff = '${mu2}'
    use_nonorthogonal_correction = false
    tke = 'TKE_p2_aux'
    rho = ${rho2}
    mu = 0.004 #'${mu}'
    alpha = 'alpha_2'
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_turb_diffusion_p2]
    type = LinearFVTurbulentMultiPhaseDiffusion
    variable = TKED_p2
    diffusion_coeff = 'mu_t_p2'
    scaling_coeff = ${sigma_eps}
    use_nonorthogonal_correction = false
    tke = 'TKE_p2_aux'
    rho = ${rho2}
    mu = ${mu2}
    alpha = 'alpha_2'
    Re_y_star = 60.0
    wall_distance = 'd'
  []
  [TKED_source_sink_p2]
    type = LinearFVTKEDMultiPhaseSourceSink
    variable = TKED_p2
    u = vel_x_p2
    v = vel_y_p2
    w = vel_z_p2
    tke = TKE_p2
    rho = ${rho2}
    mu = ${mu2}
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    wall_distance = 'd'
    alpha = 'alpha_2'
  []

  [alpha_1_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = alpha_1
    rho = 1.0 #${rho1}
    alpha = alpha_1
    MULES_iterations = ${MULES_iterations}
  []
  [alpha_1_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_1
    rhie_chow_user_object = 'rc_p1'
    c_alpha = ${c_alpha}
    rho = 1.0 #${rho1}
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    MULES_iterations = ${MULES_iterations}
  []
  [alpha_1_diff]
    type = LinearFVDiffusion
    variable = alpha_1
    diffusion_coeff = 0.1
    use_nonorthogonal_correction = false
  []

  [alpha_2_time]
    type = LinearFVMultiPhaseTimeDerivative
    variable = alpha_2
    rho = 1.0 #${rho2}
    alpha = alpha_2
    MULES_iterations = ${MULES_iterations}
  []
  [alpha_2_advection]
    type = LinearFVMultiPhaseFractionAdvection
    variable = alpha_2
    rhie_chow_user_object = 'rc_p2'
    c_alpha = ${c_alpha}
    rho = 1.0 #${rho2}
    advected_interp_method = ${advected_interp_method}
    limiter_method = ${limiter_method}
    use_nonorthogonal_correction = false
    MULES_iterations = ${MULES_iterations}
  []
  [alpha_2_diff]
    type = LinearFVDiffusion
    variable = alpha_2
    diffusion_coeff = 0.1
    use_nonorthogonal_correction = false
  []
[]

[LinearFVBCs]

  [walls-u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_x_p1
    functor = 0.0
  []
  [walls-v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_y_p1
    functor = 0.0
  []
  [walls-w_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_z_p1
    functor = 0.0
  []

  [walls-u-impeller_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_x_p1
    functor = 'vel_x_impeller'
  []
  [walls-v-impeller_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_y_p1
    functor = 'vel_y_impeller'
  []
  [walls-w-impeller_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_z_p1
    functor = 0.0
  []

  [walls-u_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_x_p2
    functor = 'vel_x_impeller'
  []
  [walls-v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_y_p2
    functor = 'vel_y_impeller'
  []
  [walls-w_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = '${walls_fixed}'
    variable = vel_z_p2
    functor = 0.0
  []

  [walls-u-impeller_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_x_p2
    functor = 'vel_x_impeller'
  []
  [walls-v-impeller_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_y_p2
    functor = 'vel_y_impeller'
  []
  [walls-w-impeller_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'walls-impeller'
    variable = vel_z_p2
    functor = 0.0
  []
  
  [walls_p]
    type = LinearFVExtrapolatedPressureBC
    boundary = '${walls_all}'
    variable = pressure
  []
  [outlet_p]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'outlet'
    variable = pressure
    functor = 1.4
  []

  [outlet_u_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x_p1
    use_two_term_expansion = false
    boundary = 'outlet'
  []
  [outlet_v_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y_p1
    use_two_term_expansion = false
    boundary = 'outlet'
  []
  [outlet_w_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_z_p1
    use_two_term_expansion = false
    boundary = 'outlet'
  []
  [inlet_u_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x_p1
    functor = 0.0
    boundary = 'inlet'
  []
  [inlet_v_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y_p1
    functor = 0.0
    boundary = 'inlet'
  []
  [inlet_w_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_z_p1
    functor = 1.0
    boundary = 'inlet'
  []

  [outlet_u_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_x_p2
    use_two_term_expansion = false
    boundary = 'outlet'
  []
  [outlet_v_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_y_p2
    use_two_term_expansion = false
    boundary = 'outlet'
  []
  [outlet_w_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    variable = vel_z_p2
    use_two_term_expansion = false
    boundary = 'outlet'
  []
  [inlet_u_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_x_p2
    functor = 0.0
    boundary = 'inlet'
  []
  [inlet_v_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_y_p2
    functor = 0.0
    boundary = 'inlet'
  []
  [inlet_w_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    variable = vel_z_p2
    functor = 1.0
    boundary = 'inlet'
  []

  [outlet_TKE_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKE_p1
    use_two_term_expansion = false
  []
  [outlet_TKED_p1]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = TKED_p1
    use_two_term_expansion = false
  []
  [inlet_TKE_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKE_p1
    functor = 0.1
  []
  [inlet_TKED_p1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKED_p1
    functor = 0.1
  []

  [outlet_TKE_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'inlet outlet'
    variable = TKE_p2
    use_two_term_expansion = false
  []
  [outlet_TKED_p2]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'inlet outlet'
    variable = TKED_p2
    use_two_term_expansion = false
  []
  [inlet_TKE_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKE_p2
    functor = 0.1
  []
  [inlet_TKED_p2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = TKED_p2
    functor = 0.1
  []

  [inlet_alpha_1]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = alpha_1
    functor = ${phase_1_inlet}
  []
  [outlet_alpha_1]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = alpha_1
    use_two_term_expansion = false
  []

  [inlet_alpha_2]
    type = LinearFVAdvectionDiffusionFunctorDirichletBC
    boundary = 'inlet'
    variable = alpha_2
    functor = ${phase_2_inlet}
  []
  [outlet_alpha_2]
    type = LinearFVAdvectionDiffusionOutflowBC
    boundary = 'outlet'
    variable = alpha_2
    use_two_term_expansion = false
  []
[]

[AuxVariables]
  [mu_t_p1]
    type = MooseLinearVariableFVReal
    initial_condition = '0.01'
  []
  [mu_t_p2]
    type = MooseLinearVariableFVReal
    initial_condition = '0.01'
  []
  [d]
    type = MooseLinearVariableFVReal
    initial_condition = '1.0'
  []
  [TKE_p1_aux]
    type = MooseLinearVariableFVReal
    initial_condition = '0.1'
  []
  [TKE_p2_aux]
    type = MooseLinearVariableFVReal
    initial_condition = '0.1'
  []
  # [alpha_2]
  #   type = MooseLinearVariableFVReal
  # []
[]

[AuxKernels]
  [compute_mu_t_p1]
    type = kEpsilonViscosityAux
    variable = mu_t_p1
    C_mu = ${C_mu}
    tke = TKE_p1
    epsilon = TKED_p1
    mu = ${mu1}
    rho = ${rho1}
    u = vel_x_p1
    v = vel_y_p1
    w = vel_z_p1
    walls = '${walls_all}'
    wall_treatment = ${wall_treatment}
    wall_distance = 'd'
    execute_on = 'NONLINEAR'
  []
  [compute_mu_t_p2]
    type = kEpsilonViscosityAux
    variable = mu_t_p2
    C_mu = ${C_mu}
    tke = TKE_p2
    epsilon = TKED_p2
    mu = ${mu2}
    rho = ${rho2}
    u = vel_x_p2
    v = vel_y_p2
    w = vel_z_p2
    walls = '${walls_all}'
    wall_treatment = ${wall_treatment}
    wall_distance = 'd'
    execute_on = 'NONLINEAR'
  []
  [compute_wall_distance]
    type = WallDistanceAux
    variable = d
    walls = '${walls_all}'
    execute_on = 'INITIAL'
  []
  [compute_TKE_p1_aux]
    type = FunctorAux
    variable = 'TKE_p1_aux'
    functor = 'TKE_p1'
    execute_on = 'NONLINEAR'
  []
  [compute_TKE_p2_aux]
    type = FunctorAux
    variable = 'TKE_p2_aux'
    functor = 'TKE_p2'
    execute_on = 'NONLINEAR'
  []
  # [populate_alpha_2]
  #   type = ParsedAux
  #   variable = alpha_2
  #   coupled_variables = 'alpha_1'
  #   expression = 'min(max(1.0 - alpha_1, 0), 1)'
  #   execute_on = 'NONLINEAR'
  # []
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

  number_of_phases = 2

  momentum_l_abs_tol = 1e-10
  pressure_l_abs_tol = 1e-10
  turbulence_l_abs_tol = 1e-14
  phase_l_abs_tol = 1e-14

  momentum_l_tol = 0
  pressure_l_tol = 0
  turbulence_l_tol = 0
  phase_l_tol = 0

  rhie_chow_user_objects = 'rc_p1 rc_p2'
  momentum_systems = 'u_system_p1 v_system_p1 w_system_p1; u_system_p2 v_system_p2 w_system_p2'
  pressure_system = 'pressure_system'
  turbulence_systems = 'TKE_system_p1 TKED_system_p1; TKE_system_p2 TKED_system_p2'
  phase_systems = 'alpha_system_p1 alpha_system_p2'

  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  turbulence_equation_relaxation = '0.25 0.25'
  turbulence_field_relaxation = '0.25 0.25'
  phase_equation_relaxation = 0.1

  num_iterations = 5
  dt = 0.01
  num_steps = 5000
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

  # Interface tratment
  enforce_phase_sum = true
  activate_interface_shapening = false
  shapening_type = 'heaviside'
  smoothing_constant = 100.0
  MULES_iterations = ${MULES_iterations}

[]

[Outputs]
  exodus = true
[]
