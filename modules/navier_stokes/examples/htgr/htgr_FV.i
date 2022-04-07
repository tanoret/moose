# ==============================================================================
# Model description
# ------------------------------------------------------------------------------
# Idaho Falls, INL, August 10, 2021 12:56 PM
# Author(s): David A. Reger, Dr. Paolo Balestra
# Idaho Falls, INL, March 7, 2022
# Author(s): Dr. Mauricio Tano
# ==============================================================================
# - ss: Steady state pseudo transient.
# ==============================================================================
# - The Model has been built based on [1].
# ------------------------------------------------------------------------------
# [1]
# ==============================================================================
#
# ==============================================================================
# MODEL PARAMETERS
# ==============================================================================
# Blocks -----------------------------------------------------------------------
riser_fluid_blocks = ' 1 2 3 '
fluid_only_blocks = ' 4 '
pbed_blocks = ' 5 '
bypass_fluid_blocks = ' 9 '
outlet_fluid_blocks = ' 6 7 8 10 '
solid_only_blocks = ' 11 12 13 14 15 16 17 '
fluid_porous_blocks = '${riser_fluid_blocks} ${pbed_blocks} ${bypass_fluid_blocks} ${outlet_fluid_blocks}'
fluid_blocks = '${fluid_porous_blocks} ${fluid_only_blocks}'

# Geometry ---------------------------------------------------------------------
pebble_bed_porosity = 0.39
void_porosity = 0.95
riser_porosity = 0.22
bypass_porosity = 0.10
outlet_channels_porosity = 0.25
pebble_bed_radius = 1.2
pebble_bed_height = 8.93 # No cone
pebble_bed_baf = 2.245 # No Cone
pebble_diameter = 0.06
inlet_free_flow_area = 4.2976

# BCs --------------------------------------------------------------------------
reactor_inlet_T_fluid = 533.25
reactor_outlet_pressure = 5.84e6
reactor_total_mfr = 78.6
power_density = 4.846688343e6
heat_capacity_multiplier = 1e-2
inlet_superficial_rho_vel = ${fparse -1* reactor_total_mfr/inlet_free_flow_area}

# ==============================================================================
# GEOMETRY AND MESH
# ==============================================================================
[Mesh]
  type = MeshGeneratorMesh

  uniform_refine = 0

  [cartesian_mesh]
    type = CartesianMeshGenerator
    # --------------------------------------------------------------------------
    # 01 Inlet
    # 02 Riser
    # 03 Diffuser
    # 04 Upper plenum
    # 05 Pebble Bed
    # 06 Bottom reflector
    # 07 Bottom plenum
    # 08 Outlet channel
    # 09 CRs He Gap
    # 10 Outlet
    # 11 12 13 14 15 16 17 Plain Graphite
    # --------------------------------------------------------------------------
    dim = 2

    dx = '0.360 0.210 0.210 0.210 0.210
          0.060
          0.100 0.360
          0.180 0.210 '

    ix = '3 2 2 2 2
          1
          1 3
          2 2'

    dy = '0.2150 0.1500 0.9500 0.3900 0.1350 0.1350 0.1350 0.1350 0.3600
          0.7791 0.7791 0.7791 0.7791 0.7791 0.7791 0.7791 0.7791 0.7791 0.7791 0.7791
          0.5500 0.8500'

    # iy = '2 2 2 2 2 2 2 2 2
    #       2 2 2 2 2 2 2 2 2 2 2
    #       2 2'
    iy = '2 2 10 4 1 1 1 1 4
          8 8 8 8 8 8 8 8 8 8 8
          5 9'

    subdomain_id = '
                    16 16 16 16 16  16  16   16  16  17

                     8  8  8  8  8   8   8    8  10  17

                    11  7  7  7  7  12   9   13  14  17

                    11  6  6  6  6  12   9   13   1  17
                     5  6  6  6  6  12   9   13   2  17
                     5  5  6  6  6  12   9   13   2  17
                     5  5  5  6  6  12   9   13   2  17
                     5  5  5  5  6  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17
                     5  5  5  5  5  12   9   13   2  17

                     4  4  4  4  4   3   3    3   3  17

                    15 15 15 15 15  15  15   15  15  17
                    '
  []
  [reactor_inlet]
    type = SideSetsAroundSubdomainGenerator
    input = cartesian_mesh
    fixed_normal = true
    normal = '1 0 0'
    new_boundary = reactor_inlet
    block = 1
  []
  [reactor_outlet]
    type = SideSetsAroundSubdomainGenerator
    input = reactor_inlet
    fixed_normal = true
    normal = '1 0 0'
    new_boundary = reactor_outlet
    block = 10
  []

  # Vertical Walls.
  [vertical_walls_01]
    type = SideSetsBetweenSubdomainsGenerator
    input = reactor_outlet
    primary_block = ' 2 3
                    '
    paired_block = ' 17

                   '
    new_boundary = vertical_walls
  []
  [vertical_walls_02]
    type = SideSetsBetweenSubdomainsGenerator
    input = vertical_walls_01
    primary_block = '
                      1 2 9
                    '
    paired_block = '
                     12 13
                   '
    new_boundary = vertical_walls
  []
  [vertical_walls_03]
    type = SideSetsBetweenSubdomainsGenerator
    input = vertical_walls_02
    primary_block = '
                      5 6 7
                    '
    paired_block = '
                     12
                   '
    new_boundary = vertical_walls
  []
  [vertical_walls_04]
    type = SideSetsBetweenSubdomainsGenerator
    input = vertical_walls_03
    primary_block = '
                      6 7
                    '
    paired_block = '
                     11
                   '
    new_boundary = vertical_walls
  []

  # Horizontal Walls.
  [horizontal_walls_01]
    type = SideSetsBetweenSubdomainsGenerator
    input = vertical_walls_04
    primary_block = '
                      3 4
                    '
    paired_block = '
                     12 13 15
                   '
    new_boundary = horizontal_walls
  []
  [horizontal_walls_02]
    type = SideSetsBetweenSubdomainsGenerator
    input = horizontal_walls_01
    primary_block = ' 8 10
                    '
    paired_block = ' 11 16 12 13 14
                   '
    new_boundary = horizontal_walls
  []
  [horizontal_walls_03]
    type = SideSetsBetweenSubdomainsGenerator
    input = horizontal_walls_02
    primary_block = ' 1 5
                    '
    paired_block = ' 14 11
                   '
    new_boundary = horizontal_walls
  []
  [side_reflector_barrel_gap]
    type = SideSetsAroundSubdomainGenerator
    input = horizontal_walls_03
    fixed_normal = true
    normal = '1 0 0'
    new_boundary = side_reflector_barrel_gap
    block = 17
  []
[]

# FV Parameters ----------------------------------------------------------------
advected_interp_method='average'
velocity_interp_method='rc'

[GlobalParams]
  fp = helium_obj
  gravity = '0 -9.81 0'
  pebble_diameter = ${pebble_diameter}
  rhie_chow_user_object = 'rc'
  advected_interp_method = ${advected_interp_method}
  velocity_interp_method = ${velocity_interp_method}
  porosity = porosity
  Darcy_name = 'Darcy_coefficient'
  Forchheimer_name = 'Forchheimer_coefficient'
  pressure = pressure
  u = v_x
  v = v_y
  rho = 1000 #'rho'
  mu = 0.1 #'mu'
  k = 0.1
  cp = 1000.0
  alpha_name = 0.001 #'beta'
  solid_rho = 1000.0
  ref_temperature = ${reactor_inlet_T_fluid}
  T_fluid = 'T_fluid'
  smoothing_layers = 2
  consistent_scaling = 0.0
[]

[Problem]
  coord_type = RZ
  rz_coord_axis = Y
  kernel_coverage_check = false
  material_coverage_check = false
[]

# ==============================================================================
# VARIABLES AND KERNELS
# ==============================================================================
[UserObjects]
  [rc]
    type = PINSFVRhieChowInterpolator
    block = '${fluid_blocks}'
  []
[]

[Variables]
  [pressure]
    type = INSFVPressureVariable
    initial_condition = ${reactor_outlet_pressure}
    scaling = 1.0
    block = '${fluid_blocks}'
  []
  [v_x]
    type = PINSFVSuperficialVelocityVariable
    scaling = 1.0
    initial_condition = 0.0
    block ='${fluid_blocks}'
  []
  [v_y]
    type = PINSFVSuperficialVelocityVariable
    scaling = 1.0
    initial_condition = 0.0
    block = '${fluid_blocks}'
  []
  [T_fluid]
    type = INSFVEnergyVariable
    scaling = 1.0e-3
    initial_condition = ${reactor_inlet_T_fluid}
    block = '${fluid_blocks}'
  []
  [T_solid]
    type = MooseVariableFVReal
    scaling = 1.0e-3
    initial_condition = ${reactor_inlet_T_fluid}
    block = '${fluid_porous_blocks} ${solid_only_blocks}'
  []
[]

[FVKernels]
  # Mass Equation.
  [mass]
    type = PINSFVMassAdvection
    variable = pressure
  []

  # Momentum x component equation.
  [vel_x_time]
    type = PINSFVMomentumTimeDerivative
    variable = v_x
    momentum_component = 'x'
  []
  [vel_x_advection]
    type = PINSFVMomentumAdvection
    variable = v_x
    momentum_component = 'x'
  []
  [vel_x_viscosity]
    type = PINSFVMomentumDiffusion
    variable = v_x
    momentum_component = 'x'
  []
  [u_pressure]
    type = PINSFVMomentumPressure
    variable = v_x
    momentum_component = 'x'
  []
  [u_friction]
    type = PINSFVMomentumFriction
    variable = v_x
    momentum_component = 'x'
  []
  [u_gravity]
    type = PINSFVMomentumGravity
    variable = v_x
    momentum_component = 'x'
  []
  [u_boussinesq]
    type = PINSFVMomentumBoussinesq
    variable = v_x
    momentum_component = 'x'
  []

  # Momentum y component equation.
  [vel_y_time]
    type = PINSFVMomentumTimeDerivative
    variable = v_y
    momentum_component = 'y'
  []
  [vel_y_advection]
    type = PINSFVMomentumAdvection
    variable = v_y
    momentum_component = 'y'
  []
  [vel_y_viscosity]
    type = PINSFVMomentumDiffusion
    variable = v_y
    momentum_component = 'y'
  []
  [v_pressure]
    type = PINSFVMomentumPressure
    variable = v_y
    momentum_component = 'y'
  []
  [v_friction]
    type = PINSFVMomentumFriction
    variable = v_y
    momentum_component = 'y'
  []
  [gravity]
    type = PINSFVMomentumGravity
    variable = v_y
    momentum_component = 'y'
  []
  [buoyancy_boussinesq]
    type = PINSFVMomentumBoussinesq
    variable = v_y
    momentum_component = 'y'
  []

#   # Fluid Energy equation.
#   [temp_fluid_time]
#     type = PINSFVEnergyTimeDerivative
#     variable = T_fluid
#     is_solid = false
#   []
#   [temp_fluid_advection]
#     type = PINSFVEnergyAdvection
#     variable = T_fluid
#     advected_quantity = 'rho_cp_temp'
#   []
#   [temp_fluid_conduction]
#     type = PINSFVEnergyDiffusion
#     variable = T_fluid
#     effective_diffusivity = false
#   []
#   [temp_solid_to_fluid]
#     type = PINSFVEnergyAmbientConvection
#     variable = T_fluid
#     is_solid = false
#     h_solid_fluid = 0.1 #'alpha'
#   []
#
#   # Solid Energy equation.
#   [temp_solid_time]
#     type = PINSFVEnergyTimeDerivative
#     variable = T_solid
#     cp = 1000.0 #'cp_s'
#     rho = ${solid_rho}
#     is_solid = true
#     block = '${fluid_porous_blocks}'
#   []
#   [temp_solid_conduction]
#     type = FVDiffusion
#     variable = T_solid
#     coeff = 10.0 #'kappa_s'
#     block = '${fluid_porous_blocks}'
#   []
#   [temp_solid_source]
#     type = FVCoupledForce
#     variable = T_solid
#     v = power_distribution
#     block = '${fluid_porous_blocks}'
#   []
#   [temp_fluid_to_solid]
#     type = PINSFVEnergyAmbientConvection
#     variable = T_solid
#     is_solid = true
#     h_solid_fluid = 0.1 #'alpha'
#     block = '${fluid_porous_blocks}'
#   []
# []

# ==============================================================================
# AUXVARIABLES AND AUXKERNELS
# ==============================================================================
[AuxVariables]
  [porosity]
    type = MooseVariableFVReal
    block = '${fluid_blocks}'
  []
  # [power_distribution]
  #   type = MooseVariableFVReal
  #   block = '${fluid_porous_blocks}'
  # []
[]

# ==============================================================================
# INITIAL CONDITIONS AND FUNCTIONS
# ==============================================================================
[ICs]
  # Porosity.
  [pbed_porosity_initialization]
    type = FunctionIC
    variable = porosity
    function = ${pebble_bed_porosity}
    block = '${pbed_blocks}'
  []
  [riser_porosity_initialization]
    type = FunctionIC
    variable = porosity
    function = ${riser_porosity}
    block = '${riser_fluid_blocks}'
  []
  [outlet_channels_porosity_initialization]
    type = FunctionIC
    variable = porosity
    function = ${outlet_channels_porosity}
    block = '${outlet_fluid_blocks}'
  []
  [cavity_porosity_initialization]
    type = FunctionIC
    variable = porosity
    function = ${void_porosity}
    block = '${fluid_only_blocks}'
  []
  [bypass_porosity_initialization]
    type = FunctionIC
    variable = porosity
    function = ${bypass_porosity}
    block = '${bypass_fluid_blocks}'
  []

  # [pow_init1]
  #   type = FunctionIC
  #   variable = power_distribution
  #   function = ${power_density}
  #   block = '${pbed_blocks}'
  # []
[]

[Functions]
[]

# ==============================================================================
# FLUID PROPERTIES, MATERIALS AND USER OBJECTS
# ==============================================================================
[FluidProperties]
  [fp]
    type = HeliumFluidProperties
  []
[]
