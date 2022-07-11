################################################################################
## Molten Salt Fast Reactor - Euratom EVOL + Rosatom MARS Design              ##
## Pronghorn Sub-Application input file                                       ##
## Relaxation towards Steady state 3D thermal hydraulics model                ##
################################################################################
# The flow in this simulation should be initialized with a previous flow
# solution (isothermal, heated or multiphysics) OR using a viscosity rampdown
# An isothermal solution can be generated using 'run_ns_initial.i', and is
# saved in 'restart/run_ns_initial_restart.e'
# A heated solution, with a flat power distribution, can be generated with
# this script 'run_ns.i', and is saved in 'restart/run_ns_restart.e'
# A coupled neutronics-coarse TH solution can be generated with
# 'run_neutronics.i', saved in 'restart/run_neutronics_ns_restart.e'
# Material properties
rho = 4284  # density [kg / m^3]  (@1000K)
cp = 1594  # specific heat capacity [J / kg / K]
drho_dT = 0.882  # derivative of density w.r.t. temperature [kg / m^3 / K]
mu = 0.0166 # viscosity [Pa s]
k = 1.7 # thermal conductivity [W / m / K]
# https://www.researchgate.net/publication/337161399_Development_of_a_control-\
# oriented_power_plant_simulator_for_the_molten_salt_fast_reactor/fulltext/5dc9\
# 5c9da6fdcc57503eec39/Development-of-a-control-oriented-power-plant-simulator-\
# for-the-molten-salt-fast-reactor.pdf
von_karman_const = 0.41
# Turbulent properties
Pr_t = 0.9 # turbulent Prandtl number
# Derived material properties
alpha = ${fparse drho_dT / rho}  # thermal expansion coefficient
# Operating parameters
T_HX = 873.15 # heat exchanger temperature [K]
# Mass flow rate tuning, for heat exchanger pressure and temperature drop
friction = 4e3  # [kg / m^4]
pump_force = -20000. # [N / m^3]
[GlobalParams]
  rhie_chow_user_object = 'ins_rhie_chow_interpolator'
[]
################################################################################
# GEOMETRY
################################################################################
[Mesh]
  [restart]
    type = FileMeshGenerator
    use_for_exodus_restart = true
    # Depending on the file chosen, the initialization of variables should be
    # adjusted. The following variables can be initalized:
    # - vel_x, vel_y, p from isothermal simulation
    # file = 'restart/run_ns_initial_restart.e'
    # Below are initialization points created from this input file
    # The variable IC should be set from_file_var for temperature and precursors
    # - vel_x, vel_y, p, T_fluid, c_i from cosine heated simulation
    # file = 'restart/run_ns_restart.e'
    # - vel_x, vel_y, p, T_fluid, c_i from coupled multiphysics simulation
    file = '../restart/run_ns_coupled_restart.e'
  []
  [hx_top]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y > 0'
    included_subdomains = '3'
    included_neighbors = '1'
    fixed_normal = true
    normal = '0 1 0'
    new_sideset_name = 'hx_top'
    input = 'restart'
  []
  [hx_bot]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'y <-0.6'
    included_subdomains = '3'
    included_neighbors = '1'
    fixed_normal = true
    normal = '0 -1 0'
    new_sideset_name = 'hx_bot'
    input = 'hx_top'
  []
[]
[Problem]
  coord_type = 'RZ'
  rz_coord_axis = Y
[]
################################################################################
# EQUATIONS: VARIABLES, KERNELS, BOUNDARY CONDITIONS
################################################################################
[Modules]
  [NavierStokesFV]
    # General parameters
    compressibility = 'incompressible'
    add_energy_equation = true
    boussinesq_approximation = true
    # Variables, defined below for the Exodus restart
    velocity_variable = 'vel_x vel_y'
    pressure_variable = 'pressure'
    fluid_temperature_variable = 'T_fluid'
    # Material properties
    density = ${rho}
    dynamic_viscosity = ${mu}
    thermal_conductivity = ${k}
    specific_heat = 'cp'
    thermal_expansion = ${alpha}
    # Boussinesq parameters
    gravity = '0 -9.81 0'
    ref_temperature = ${T_HX}
    # Boundary conditions
    wall_boundaries = 'shield_wall reflector_wall fluid_symmetry'
    momentum_wall_types = 'wallfunction wallfunction symmetry'
    energy_wall_types = 'heatflux heatflux heatflux'
    energy_wall_function = '0 0 0'
    # Pressure pin for incompressible flow
    pin_pressure = true
    pinned_pressure_type = average
    pinned_pressure_value = 1e5
    # Turbulence parameters
    turbulence_handling = 'mixing-length'
    turbulent_prandtl = ${Pr_t}
    von_karman_const = ${von_karman_const}
    mixing_length_delta = 0.1
    mixing_length_walls = 'shield_wall reflector_wall'
    mixing_length_aux_execute_on = 'initial'
    # Numerical scheme
    momentum_advection_interpolation = 'upwind'
    mass_advection_interpolation = 'upwind'
    energy_advection_interpolation = 'upwind'
    # Heat source
    external_heat_source = power_density
    # Heat exchanger
    friction_blocks = 'hx'
    friction_types = 'FORCHHEIMER'
    friction_coeffs = ${friction}
    ambient_convection_blocks = 'hx'
    ambient_convection_alpha = ${fparse 600 * 20e3} # HX specifications
    ambient_temperature = ${T_HX}
  []
[]
[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    block = 'fuel pump hx'
    initial_from_file_var = vel_x
  []
  [vel_y]
    type = INSFVVelocityVariable
    block = 'fuel pump hx'
    initial_from_file_var = vel_y
  []
  [pressure]
    type = INSFVPressureVariable
    block = 'fuel pump hx'
    initial_from_file_var = pressure
  []
  [T_fluid]
    type = INSFVEnergyVariable
    block = 'fuel pump hx'
    initial_condition = ${T_HX}
    # initial_from_file_var = T_fluid
  []
[]
[AuxVariables]
  [power_density]
    type = MooseVariableFVReal
    block = 'fuel pump hx'
    # Power density is re-initalized by a transfer from neutronics
    [InitialCondition]
      type = FunctionIC
      function = 'cosine_guess'
      scaling_factor = ${fparse 3e9/2.81543}
    []
  []
  [ax_var]
    type = MooseVariableFVReal
  []
  [ay_var]
    type = MooseVariableFVReal
  []
[]
[AuxKernels]
  [ax_fill]
    type = ADFunctorElementalAux
    variable = ax_var
    functor = 'ax'
    execute_on = 'TIMESTEP_END'
  []
  [ay_fill]
    type = ADFunctorElementalAux
    variable = ay_var
    functor = 'ay'
    execute_on = 'TIMESTEP_END'
  []
[]
[Functions]
  # Guess to have a 3D power distribution
  [cosine_guess]
    type = ParsedFunction
    value = 'max(0, cos(x*pi/2/1.2))*max(0, cos(y*pi/2/1.1))'
  []
[]
[FVKernels]
  [pump]
    type = INSFVBodyForce
    variable = vel_y
    functor = ${pump_force}
    block = 'pump'
    momentum_component = 'y'
  []
[]
################################################################################
# MATERIALS
################################################################################
[Materials]
  # Most of these constants could be specified directly to the action
  [mu]
    type = ADGenericFunctorMaterial
    prop_names = 'mu'
    prop_values = '${mu}'
    block = 'fuel pump hx'
  []
  # [not_used]
  #   type = ADGenericFunctorMaterial
  #   prop_names = 'not_used'
  #   prop_values = 0
  #   block = 'shield reflector'
  # []
  [cp]
    type = ADGenericFunctorMaterial
    prop_names = 'cp dcp_dt'
    prop_values = '${cp} 0'
    block = 'fuel pump hx'
  []
[]
################################################################################
# EXECUTION / SOLVE
################################################################################
[Functions]
  [dts]
    type = PiecewiseConstant
    x = '0    100'
    y = '0.75 2.5'
  []
[]
[Executioner]
  type = Transient
  # Time stepping parameters
  start_time = 0.0
  end_time = 200
  # end_time will depend on the restart file chosen
  # though steady state detection can also be used
  # from _initial/no heating : 150 - 200s enough
  # from _ns/_ns_coupled/heated: 10s enough
  [TimeStepper]
    # This time stepper makes the time step grow exponentially
    # It can only be used with proper initialization
    type = IterationAdaptiveDT
    dt = 1  # chosen to obtain convergence with first coupled iteration
    growth_factor = 2
  []
  # [TimeStepper]
  #   type = FunctionDT
  #   function = dts
  # []
  steady_state_detection  = true
  steady_state_tolerance  = 1e-8
  steady_state_start_time = 10
  # Solver parameters
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -ksp_gmres_restart'
  petsc_options_value = 'lu NONZERO 50'
  line_search = 'none'
  nl_rel_tol = 1e-9
  nl_abs_tol = 2e-8
  nl_max_its = 15
  l_max_its = 50
  automatic_scaling = true
  # resid_vs_jac_scaling_param = 1
[]
################################################################################
# SIMULATION OUTPUTS
################################################################################
[Outputs]
  csv = true
  exodus = true
[]
[Postprocessors]
[]
