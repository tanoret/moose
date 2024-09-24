steel_thickness = 0.00254
pipe_D = 0.0254
width = ${fparse 10.0*pipe_D}
height = ${fparse 20.0*pipe_D}

p_outlet = 1e5
T_Salt_initial = 900

heat_source = 1e6
hx_alpha = 1e5

cp_steel = 500.0 # (J/(kg.K)) specific heat of steel
rho_steel = 8000.0 # (kg/(m3)) density of steel
k_steel = 15.0 # # (W/(m.k)) density of steel
bulk_hx = 100.0

# F = 96485.3321
# epsilon_0 = 8.854e-12
# epsilon_r = 10.0
# Foeps = ${fparse F/(epsilon_0*epsilon_r)}
Foeps = 1e6

pipe_cells = 10
elbow_cells = 25
width_cells = ${fparse 4 * elbow_cells}
height_cells = ${fparse 8 * elbow_cells}

# Steel composition
c_Cr_steel = 0.17
c_Ni_steel = 0.12
c_Fe_steel = ${fparse 1.0 - c_Cr_steel - c_Ni_steel}

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${steel_thickness} ${pipe_D} ${steel_thickness} ${width} ${steel_thickness} ${pipe_D} ${steel_thickness}'
    dy = '${steel_thickness} ${pipe_D} ${steel_thickness} ${height} ${steel_thickness} ${pipe_D} ${steel_thickness}'
    ix = '${pipe_cells} ${elbow_cells} ${pipe_cells} ${width_cells} ${pipe_cells} ${elbow_cells} ${pipe_cells}'
    iy = '${pipe_cells} ${elbow_cells} ${pipe_cells} ${height_cells} ${pipe_cells} ${elbow_cells} ${pipe_cells}'
    subdomain_id = '1 1 1 1 1 1 1
                    1 2 2 2 2 2 1
                    1 3 1 1 1 4 1
                    1 3 1 5 1 4 1
                    1 3 1 1 1 4 1
                    1 2 2 2 2 2 1
                    1 1 1 1 1 1 1'
  []
  [solid_liquid_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = cmg
    new_boundary = 'solid_liquid_interface'
    primary_block = '1'
    paired_block = '2 3 4'
  []
  [delte_internal_block]
    type = BlockDeletionGenerator
    input = solid_liquid_interface
    block = '5'
    new_boundary = 'wall_in'
  []
  [create_wall]
    type = RenameBoundaryGenerator
    input = delte_internal_block
    old_boundary = 'wall_in left top bottom right'
    new_boundary = 'wall wall wall wall wall'
  []
  [rename_blocks]
    type = RenameBlockGenerator
    input = create_wall
    old_block = '1 2 3 4'
    new_block = 'pipe flow flow-heated flow-cooled'
  []
[]

fluid_blocks = 'flow flow-heated flow-cooled'
solid_blocks = 'pipe'

[Problem]
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Variables]
  [vel_x]
    type = INSFVVelocityVariable
    initial_condition = 1e-8
    block = ${fluid_blocks}
  []
  [vel_y]
    type = INSFVVelocityVariable
    initial_condition = 1e-8
    block = ${fluid_blocks}
  []
  [pressure]
    type = INSFVPressureVariable
    initial_condition = ${p_outlet}
    block = ${fluid_blocks}
  []
  [T_fluid]
    type = INSFVEnergyVariable
    initial_condition = ${T_Salt_initial}
    block = ${fluid_blocks}
  []
  [T_solid]
    type = INSFVEnergyVariable
    initial_condition = ${T_Salt_initial}
    block = ${solid_blocks}
  []

  [Fe_s]
    type = MooseVariableFVReal
    initial_condition = ${c_Fe_steel}
    block = ${solid_blocks}
  []
  [Fe]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = ${fluid_blocks}
  []
  [Cr_s]
    type = MooseVariableFVReal
    initial_condition = ${c_Cr_steel}
    block = ${solid_blocks}
  []
  [Cr]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = ${fluid_blocks}
  []
  [Ni_s]
    type = MooseVariableFVReal
    initial_condition = ${c_Ni_steel}
    block = ${solid_blocks}
  []
  [Ni]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = ${fluid_blocks}
  []
  [phi]
    type = MooseVariableFVReal
    initial_condition = 3.0
  []
[]

[FluidProperties]
  [fluid_properties_obj]
    type = FlibeFluidProperties
  []
[]

[Modules]
  [NavierStokesFV]
    # Basic settings - weakly-compressible, turbulent flow with buoyancy
    block = ${fluid_blocks}
    compressibility = 'weakly-compressible'
    add_energy_equation = true
    gravity = '0.0 -9.81 0.0'

    # Variable naming
    velocity_variable = 'vel_x vel_y'
    pressure_variable = 'pressure'
    fluid_temperature_variable = 'T_fluid'

    # Numerical schemes
    pressure_face_interpolation = average
    momentum_advection_interpolation = upwind
    mass_advection_interpolation = upwind
    energy_advection_interpolation = upwind
    velocity_interpolation = rc

    # fluid properties
    density = 'rho'
    dynamic_viscosity = 'mu'
    thermal_conductivity = 'k'
    specific_heat = 'cp'

    # Boundary Conditions
    wall_boundaries = 'solid_liquid_interface'
    momentum_wall_types = 'noslip'
    energy_wall_types = 'heatflux'
    energy_wall_function = '0'

    # Constrain Pressure
    pin_pressure = true
    pinned_pressure_value = ${p_outlet}
    pinned_pressure_type = average-uo

    # Heat sources
    external_heat_source = 'power_source'
    ambient_convection_alpha = ${hx_alpha}
    ambient_convection_blocks = 'flow-cooled'
    ambient_temperature = ${T_Salt_initial}
  []
[]

[AuxVariables]
  [power_source]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = ${fluid_blocks}
  []
[]

[AuxKernels]
  [populate_power_source]
    type = ParsedAux
    variable = power_source
    expression = '${heat_source}'
    block = 'flow-heated'
  []
[]

[FVKernels]
  # Kernels for solve in the solid blocks
  [heat_time_solid]
    type = INSFVEnergyTimeDerivative
    variable = T_solid
    dh_dt = dh_dt
    rho = ${rho_steel}
  []
  [heat_diffusion_solid]
    type = FVDiffusion
    variable = T_solid
    coeff = ${k_steel}
  []

  [Fe_s_time]
    type = FVTimeKernel
    variable = 'Fe_s'
  []
  [diffusion_solid_Fe]
    type = FVDiffusion
    variable = 'Fe_s'
    coeff = 'diff_coef_species'
    block = ${solid_blocks}
  []
  [electrophoresis_Fe_s]
    type = FVElectrophoresisSource
    variable = 'Fe_s'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
  []
  [Fe_time]
    type = FVTimeKernel
    variable = 'Fe'
  []
  [diffusion_liquid_Fe]
    type = FVDiffusion
    variable = 'Fe'
    coeff = 'diff_coef_species'
    block = ${fluid_blocks}
  []
  [Fe_advection]
    type = INSFVScalarFieldAdvection
    variable = 'Fe'
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
    block = ${fluid_blocks}
  []
  [electrophoresis_Fe]
    type = FVElectrophoresisSource
    variable = 'Fe'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
  []

  [Cr_s_time]
    type = FVTimeKernel
    variable = 'Cr_s'
  []
  [diffusion_solid_Cr]
    type = FVDiffusion
    variable = 'Cr_s'
    coeff = 'diff_coef_species'
    block = ${solid_blocks}
  []
  [electrophoresis_Cr_s]
    type = FVElectrophoresisSource
    variable = 'Cr_s'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
  []
  [Cr_time]
    type = FVTimeKernel
    variable = 'Cr'
  []
  [diffusion_liquid_Cr]
    type = FVDiffusion
    variable = 'Cr'
    coeff = 'diff_coef_species'
    block = ${fluid_blocks}
  []
  [Cr_advection]
    type = INSFVScalarFieldAdvection
    variable = 'Cr'
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
    block = ${fluid_blocks}
  []
  [electrophoresis_Cr]
    type = FVElectrophoresisSource
    variable = 'Cr'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
  []

  [Ni_s_time]
    type = FVTimeKernel
    variable = 'Ni_s'
  []
  [diffusion_solid_Ni]
    type = FVDiffusion
    variable = 'Ni_s'
    coeff = 'diff_coef_species'
    block = ${solid_blocks}
  []
  [electrophoresis_Ni_s]
    type = FVElectrophoresisSource
    variable = 'Ni_s'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
  []
  [Ni_time]
    type = FVTimeKernel
    variable = 'Ni'
  []
  [diffusion_liquid_Ni]
    type = FVDiffusion
    variable = 'Ni'
    coeff = 'diff_coef_species'
    block = ${fluid_blocks}
  []
  [Ni_advection]
    type = INSFVScalarFieldAdvection
    variable = 'Ni'
    rhie_chow_user_object = 'ins_rhie_chow_interpolator'
    block = ${fluid_blocks}
  []
  [electrophoresis_Ni]
    type = FVElectrophoresisSource
    variable = 'Ni'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
  []


  [laplace_phi]
    type = FVDiffusion
    variable = 'phi'
    coeff = '1.0'
  []
  [Cr_source]
    type = FVCoupledForce
    variable = 'phi'
    v = 'Cr'
    coef = '${fparse 2.0*Foeps}'
    block = ${fluid_blocks}
  []
  [Fe_source]
    type = FVCoupledForce
    variable = 'phi'
    v = 'Fe'
    coef = '${fparse 2.0*Foeps}'
    block = ${fluid_blocks}
  []
  [Ni_source]
    type = FVCoupledForce
    variable = 'phi'
    v = 'Ni'
    coef = '${fparse 2.0*Foeps}'
    block = ${fluid_blocks}
  []
[]

[FVBCs]
  [phi]
    type = FVDirichletBC
    variable = 'phi'
    boundary = 'wall'
    value = 3.0
  []
[]

[FVInterfaceKernels]
  # Conjugated heat transfer with core barrel
  [convection]
    type = FVConvectionCorrelationInterface
    variable1 = T_fluid
    variable2 = T_solid
    boundary = 'solid_liquid_interface'
    h = ${bulk_hx}
    T_solid = T_solid
    T_fluid = T_fluid
    subdomain1 = ${fluid_blocks}
    subdomain2 = ${solid_blocks}
    wall_cell_is_bulk = true
  []
  [interface_Cr]
    type = FVButlerBolmerInterface
    variable1 = 'Cr'
    variable2 = 'Cr_s'
    c = 'Cr'
    c_solid = 'Cr_s'
    c_E = '-3.9'
    c_dE_dT = '-0.02'
    T_ref = 900.0
    c_Z = '2'
    c_k0 = 1e-14
    boundary = 'solid_liquid_interface'
    phi = 'phi'
    temperature = 'T_fluid'
    subdomain1 = ${fluid_blocks}
    subdomain2 = ${solid_blocks}
    wall_cell_is_bulk = true
  []
  [interface_Fe]
    type = FVButlerBolmerInterface
    variable1 = 'Fe'
    variable2 = 'Fe_s'
    c = 'Fe'
    c_solid = 'Fe_s'
    c_E = '-3.5'
    c_dE_dT = '-0.02'
    T_ref = 900.0
    c_Z = '2'
    c_k0 = 1e-14
    boundary = 'solid_liquid_interface'
    phi = 'phi'
    temperature = 'T_fluid'
    subdomain1 = ${fluid_blocks}
    subdomain2 = ${solid_blocks}
    wall_cell_is_bulk = true
  []
  [interface_Ni]
    type = FVButlerBolmerInterface
    variable1 = 'Ni'
    variable2 = 'Ni_s'
    c = 'Ni'
    c_solid = 'Ni_s'
    c_E = '-3.1'
    c_dE_dT = '-0.02'
    T_ref = 900.0
    c_Z = '2'
    c_k0 = 1e-14
    boundary = 'solid_liquid_interface'
    phi = 'phi'
    temperature = 'T_fluid'
    subdomain1 = ${fluid_blocks}
    subdomain2 = ${solid_blocks}
    wall_cell_is_bulk = true
  []
[]

[FunctorMaterials]
  [fluid_props_to_mat_props]
    type = GeneralFunctorFluidProps
    pressure = 'pressure'
    T_fluid = 'T_fluid'
    speed = 'speed'
    characteristic_length = 1.0
    fp = fluid_properties_obj
    porosity = 1.0
    block = ${fluid_blocks}
  []
  [dh_dt_mat]
    type = INSFVEnthalpyFunctorMaterial
    rho = ${rho_steel}
    temperature = T_solid
    cp = ${cp_steel}
    block = ${solid_blocks}
  []
  [diff_coef_species]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'diff_coef_species'
    subdomain_to_prop_value = 'pipe 1e-10
                               flow 1e-6
                               flow-heated 1e-6
                               flow-cooled 1e-6'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'
  automatic_scaling = true
  nl_rel_tol = 1e-3
  nl_abs_tol = 1e-5
  nl_max_its = 20

  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 10
    iteration_window = 2
    growth_factor = 2
    cutback_factor = 0.5
  []

  end_time = 1e10
  # steady_state_detection = true
  # steady_state_tolerance = 1e-12
[]

[Outputs]
  exodus = true
[]
