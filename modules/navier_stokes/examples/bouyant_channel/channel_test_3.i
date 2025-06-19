mu = 0.001
rho = 1.0
k = 5.0
cp = 4816
alpha = 150
alpha_b = 0.001
Pr_t = 0.9
advected_interp_method = 'average'
velocity_interp_method = 'rc'
g = 0.0

pressure_tag = "pressure_grad"

H = 1 #halfwidth of the channel
bulk_u = 1

### k-epsilon Closure Parameters ###
sigma_k = 1.0
sigma_eps = 1.3
C1_eps = 1.44
C2_eps = 1.92
C_mu = 0.09

### Initial and Boundary Conditions ###
intensity = 0.01
k_init = '${fparse 1.5*(intensity * bulk_u)^2}'
eps_init = '${fparse C_mu^0.75 * k_init^1.5 / H}'

### Modeling parameters ###
walls = 'top bottom'
wall_treatment = 'eq_newton' # Options: eq_newton, eq_incremental, eq_linearized, neq

[Mesh]
  [mesh]
    type = CartesianMeshGenerator
    dim = 2
    dx = '10.0'
    dy = '1.0'
    ix = '500'
    iy = '50'
  []
[]

[GlobalParams]
  rhie_chow_user_object = 'rc'
[]

[Problem]
  nl_sys_names = 'u_system v_system pressure_system energy_system  TKE_system TKED_system'
  previous_nl_solution_required = true
[]

[UserObjects]
  [rc]
    type = INSFVRhieChowInterpolatorSegregated
    u = vel_x
    v = vel_y
    pressure = pressure
  []
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    initial_condition = 0.5
    solver_sys = u_system
    two_term_boundary_expansion = false
  []
  [vel_y]
    type = INSFVVelocityVariable
    initial_condition = 0.0
    solver_sys = v_system
    two_term_boundary_expansion = false
  []
  [pressure]
    type = INSFVPressureVariable
    solver_sys = pressure_system
    initial_condition = 0.2
    two_term_boundary_expansion = false
  []
  [T_fluid]
    type = INSFVEnergyVariable
    initial_condition = 300
    solver_sys = energy_system
    two_term_boundary_expansion = false
  []
  [TKE]
    type = INSFVEnergyVariable
    solver_sys = TKE_system
    initial_condition = ${k_init}
  []
  [TKED]
    type = INSFVEnergyVariable
    solver_sys = TKED_system
    initial_condition = ${eps_init}
  []
[]

[FVKernels]
  [u_advection]
    type = INSFVMomentumAdvection
    variable = vel_x
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'x'
  []
  [u_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_x
    mu = ${mu}
    momentum_component = 'x'
  []
  [u_viscosity_turbulent]
    type = INSFVTurbulentAnisotropyDiffusion
    variable = vel_x
    coeff = 1
    bc0 = b_00_torch_func
    bc1 = b_01_torch_func
    rho = ${rho}
    k = TKE 
    #momentum_component = 'x'
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = vel_x
    momentum_component = 'x'
    pressure = pressure
    extra_vector_tags = ${pressure_tag}
  []
  [u_gravity]
    type = INSFVMomentumGravity
    momentum_component = 'x'
    variable = vel_x
    gravity = '0 ${g} 0'
    rho = ${rho}
  []

  [v_advection]
    type = INSFVMomentumAdvection
    variable = vel_y
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    rho = ${rho}
    momentum_component = 'y'
  []
  [v_viscosity]
    type = INSFVMomentumDiffusion
    variable = vel_y
    mu = ${mu}
    momentum_component = 'y'
  []
  [v_viscosity_turbulent]
    type = INSFVTurbulentAnisotropyDiffusion
    variable = vel_y
    coeff = 1
    bc0 = b_10_torch_func
    bc1 = b_11_torch_func
    rho = ${rho}
    k = TKE
    #momentum_component = 'y'
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = vel_y
    momentum_component = 'y'
    pressure = pressure
    extra_vector_tags = ${pressure_tag}
  []
  [v_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = vel_y
    momentum_component = 'y'
    T_fluid = T_fluid
    gravity = '0 ${g} 0'
    rho = ${rho}
    ref_temperature = 300.0
  []
  [v_gravity]
    type = INSFVMomentumGravity
    variable = vel_y
    momentum_component = 'y'
    gravity = '0 ${g} 0'
    rho = ${rho}
  []

  [p_diffusion]
    type = FVAnisotropicDiffusion
    variable = pressure
    coeff = "Ainv"
    coeff_interp_method = 'average'
  []
  [p_source]
    type = FVDivergence
    variable = pressure
    vector_field = "HbyA"
    force_boundary_execution = true
  []

  [energy_advection]
    type = INSFVEnergyAdvection
    variable = T_fluid
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
  []
  [energy_diffusion]
    type = FVDiffusion
    coeff = ${k}
    variable = T_fluid
  []
  [temp_turb_conduction]
    type = FVDiffusion
    coeff = 'k_t'
    variable = T_fluid
  []

  [TKE_advection]
    type = INSFVTurbulentAdvection
    variable = TKE
    rho = ${rho}
  []
  [TKE_diffusion]
    type = INSFVTurbulentDiffusion
    variable = TKE
    coeff = ${mu}
  []
  [TKE_diffusion_turbulent]
    type = INSFVTurbulentDiffusion
    variable = TKE
    coeff = 'mu_t_torch_func'
    scaling_coef = ${sigma_k}
    type = INSFVTurbulentAnisotropyDiffusion
    variable = TKE
    coeff = ${sigma_k}
    bc0 = b_00_torch_func
    bc1 = b_01_torch_func
    rho = ${rho}
    k = TKE 
    #momentum_component = 'x'
  []
  [TKE_source_sink]
    type = INSFVTKESourceSink
    variable = TKE
    u = vel_x
    v = vel_y
    epsilon = TKED
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t_torch_func'
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    anisotropy_corrections = true
    b_name = 'b_00_torch_func b_01_torch_func b_10_torch_func b_11_torch_func'
  []

  [TKED_advection]
    type = INSFVTurbulentAdvection
    variable = TKED
    rho = ${rho}
    walls = ${walls}
  []
  [TKED_diffusion]
    type = INSFVTurbulentDiffusion
    variable = TKED
    coeff = ${mu}
    walls = ${walls}
  []
  [TKED_diffusion_turbulent]
    type = INSFVTurbulentDiffusion
    variable = TKED
    coeff = 'mu_t_torch_func'
    scaling_coef = ${sigma_eps}
    walls = ${walls}
  []
  [TKED_source_sink]
    type = INSFVTKEDSourceSink
    variable = TKED
    u = vel_x
    v = vel_y
    k = TKE
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t_torch_func'
    C1_eps = ${C1_eps}
    C2_eps = ${C2_eps}
    walls = ${walls}
    wall_treatment = ${wall_treatment}
    anisotropy_corrections = true
    b_name = 'b_00_torch_func b_01_torch_func b_10_torch_func b_11_torch_func'
  []
[]

[FVBCs]
  [inlet-u]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_x
    function = '1.0'
  []
  [inlet-v]
    type = INSFVInletVelocityBC
    boundary = 'left'
    variable = vel_y
    function = '0.0'
  []
  [walls-u]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_x
    function = 0.0
  []
  [walls-v]
    type = INSFVNoSlipWallBC
    boundary = 'top bottom'
    variable = vel_y
    function = 0.0
  []
  [outlet_p]
    type = INSFVOutletPressureBC
    boundary = 'right'
    variable = pressure
    function = 1.0
  []
  [zero-grad-pressure]
    type = FVFunctionNeumannBC
    variable = pressure
    boundary = 'top left bottom'
    function = 0.0
  []
  [inlet_t]
    type = FVDirichletBC
    boundary = 'left'
    variable = T_fluid
    value = 300
  []
  [top_t]
    type = FVDirichletBC
    boundary = 'top'
    variable = T_fluid
    value = 400
  []
  [bottom_t]
    type = FVDirichletBC
    boundary = 'bottom'
    variable = T_fluid
    value = 300
  []
  [inlet_TKE]
    type = INSFVInletIntensityTKEBC
    boundary = 'left'
    variable = TKE
    u = vel_x
    v = vel_y
    intensity = ${intensity}
  []
  [inlet_TKED]
    type = INSFVMixingLengthTKEDBC
    boundary = 'left'
    variable = TKED
    k = TKE
    characteristic_length = '${fparse 2*H}'
  []
  [walls_mu_t]
    type = INSFVTurbulentViscosityWallFunction
    boundary = 'top bottom'
    variable = mu_t_torch_func
    u = vel_x
    v = vel_y
    rho = ${rho}
    mu = ${mu}
    mu_t = 'mu_t_torch_func'
    k = TKE
    wall_treatment = ${wall_treatment}
  []
  ##########################################################
[]

[Executioner]
  type = SIMPLENonlinearAssembly
  momentum_l_abs_tol = 1e-11
  pressure_l_abs_tol = 1e-11
  energy_l_abs_tol = 1e-11
  momentum_l_tol = 0
  pressure_l_tol = 0
  energy_l_tol = 0
  rhie_chow_user_object = 'rc'
  momentum_systems = 'u_system v_system'
  pressure_system = 'pressure_system'
  energy_system = 'energy_system'
  turbulence_systems = 'TKED_system TKE_system'
  pressure_gradient_tag = ${pressure_tag}
  momentum_equation_relaxation = 0.7
  pressure_variable_relaxation = 0.3
  energy_equation_relaxation = 0.5
  turbulence_equation_relaxation = '0.25 0.25'
  num_iterations = 1000
  pressure_absolute_tolerance = 1e-10
  momentum_absolute_tolerance = 1e-10
  energy_absolute_tolerance = 1e-10
  turbulence_absolute_tolerance = '1e-10 1e-10'
  print_fields = false
  continue_on_max_its = true
[]

[UserObjects]
  [cody_net]
    type = TorchScriptUserObject
    filename = "cody_net.pt"
    load_during_construction = true
    execute_on = INITIAL
  []
[]

[AuxVariables]
  [b_00_torch_func]
    type = MooseVariableFVReal
    initial_condition = 0.3
  []
  [b_01_torch_func]
    type = MooseVariableFVReal
    initial_condition = 0.3
  []
  [b_10_torch_func]
    type = MooseVariableFVReal
    initial_condition = 0.3
  []
  [b_11_torch_func]
    type = MooseVariableFVReal
    initial_condition = 0.3
  []
  [mu_t_torch_func_temp]
    type = MooseVariableFVReal
    initial_condition = 0.1
  []
  [mu_t_torch_func]
    type = MooseVariableFVReal
    initial_condition = 0.1
  []  
[]

[AuxKernels]
 [populate_b_00_torch_func]
   type = MaterialRealAux
   variable = b_00_torch_func
   property = 'b_00'
 []
[populate_b_01_torch_func]
   type = MaterialRealAux
   variable = b_01_torch_func
   property = 'b_01'
 []
[populate_b_10_torch_func]
   type = MaterialRealAux
   variable = b_10_torch_func
   property = 'b_10'
 []
[populate_b_11_torch_func]
   type = MaterialRealAux
   variable = b_11_torch_func
   property = 'b_11'
 []
 [populate_mu_t_torch_func_temp]
  type = MaterialRealAux
  variable = mu_t_torch_func_temp
  property = 'b_00'
 []
 [populate_mu_t_torch_func]
  type = ParsedAux
  variable = mu_t_torch_func
  coupled_variables = mu_t_torch_func_temp
  expression = '- ${rho} * mu_t_torch_func_temp'
 []
[]

[Materials]
  [const_functor]
    type = ADGenericFunctorMaterial
    prop_names = 'cp alpha alpha_b'
    prop_values = '${cp} ${alpha} ${alpha_b}'
  []
  #[const_mu_t]
  #  type = ADGenericFunctorMaterial
  #  prop_names = 'mu_t'
  #  prop_values = '0.1'
  #[]
  [net_material] # Populates anisotropy
   type = TorchScriptTurbulentAnisotropyMaterial
   torch_script_userobject = cody_net
   u = vel_x
   v = vel_y
   k = TKE
   eps = TKED
   debug = false
   use_NN = false
  []
  [k_t]
    type = ADParsedFunctorMaterial
    expression = 'mu_t_torch_func * cp / Pr_t'
    functor_names = 'mu_t_torch_func ${cp} ${Pr_t}'
    functor_symbols = 'mu_t_torch_func cp Pr_t'
    property_name = 'k_t'
  []
  [ins_fv]
    type = INSFVEnthalpyFunctorMaterial
    rho = ${rho}
    temperature = 'T_fluid'
  []
[]

[Outputs]
  exodus = true
  csv = false
  perf_graph = false
  print_nonlinear_residuals = false
  print_linear_residuals = true
[]
