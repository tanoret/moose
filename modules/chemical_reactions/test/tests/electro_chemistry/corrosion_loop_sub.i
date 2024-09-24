steel_thickness = 0.00254
pipe_D = 0.0254
width = ${fparse 10.0*pipe_D}
height = ${fparse 20.0*pipe_D}

# F = 96485.3321
# epsilon_0 = 8.854e-12
# epsilon_r = 10.0
# Foeps = ${fparse F/(epsilon_0*epsilon_r)}
Foeps = 10000.0

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '${steel_thickness} ${pipe_D} ${steel_thickness} ${width} ${steel_thickness} ${pipe_D} ${steel_thickness}'
    dy = '${steel_thickness} ${pipe_D} ${steel_thickness} ${height} ${steel_thickness} ${pipe_D} ${steel_thickness}'
    ix = '3 10 3 40 3 10 3'
    iy = '3 10 3 80 3 10 3'
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
  # kernel_coverage_check = false
  # material_coverage_check = false
[]

[Variables]
  [Fe_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = ${solid_blocks}
  []
  [Fe]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = ${fluid_blocks}
  []
  [Cr_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = ${solid_blocks}
  []
  [Cr]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = ${fluid_blocks}
  []
  [phi]
    type = MooseVariableFVReal
    initial_condition = 0.0
  []
[]

[AuxVariables]
  [T_fluid]
    type = MooseVariableFVReal
    initial_condition = 900.0
  []
[]

[FVKernels]
  [Fe_s_time]
    type = FVTimeKernel
    variable = 'Fe_s'
    block = ${solid_blocks}
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
    block = ${solid_blocks}
  []
  [Fe_time]
    type = FVTimeKernel
    variable = 'Fe'
    block = ${fluid_blocks}
  []
  [diffusion_liquid_Fe]
    type = FVDiffusion
    variable = 'Fe'
    coeff = 'diff_coef_species'
    block = ${fluid_blocks}
  []
  # [Fe_advection]
  #   type = INSFVScalarFieldAdvection
  #   variable = 'Fe'
  #   rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  #   block = ${fluid_blocks}
  # []
  [electrophoresis_Fe]
    type = FVElectrophoresisSource
    variable = 'Fe'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
    block = ${fluid_blocks}
  []

  [Cr_s_time]
    type = FVTimeKernel
    variable = 'Cr_s'
    block = ${solid_blocks}
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
    block = ${solid_blocks}
  []
  [Cr_time]
    type = FVTimeKernel
    variable = 'Cr'
    block = ${fluid_blocks}
  []
  [diffusion_liquid_Cr]
    type = FVDiffusion
    variable = 'Cr'
    coeff = 'diff_coef_species'
    block = ${fluid_blocks}
  []
  # [Cr_advection]
  #   type = INSFVScalarFieldAdvection
  #   variable = 'Cr'
  #   rhie_chow_user_object = 'ins_rhie_chow_interpolator'
  #   block = ${fluid_blocks}
  # []
  [electrophoresis_Cr]
    type = FVElectrophoresisSource
    variable = 'Cr'
    temperature = 'T_fluid'
    D = 'diff_coef_species'
    phi = 'phi'
    block = ${fluid_blocks}
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
    coef = '${fparse Foeps}'
    block = ${fluid_blocks}
  []
  [Fe_source]
    type = FVCoupledForce
    variable = 'phi'
    v = 'Fe'
    coef = '${fparse Foeps}'
    block = ${fluid_blocks}
  []
[]

[FVBCs]
  [phi]
    type = FVDirichletBC
    variable = 'phi'
    boundary = 'wall'
    value = 0.0
  []
[]

[FVInterfaceKernels]
  [interface_Cr]
    type = FVButlerBolmerInterface
    variable1 = 'Cr'
    variable2 = 'Cr_s'
    c = 'Cr'
    c_solid = 'Cr_s'
    c_E = '-0.9'
    c_dE_dT = '-0.0005'
    c_Z = '2'
    c_k0 = 1e-15
    boundary = 'solid_liquid_interface'
    phi = 0.0 #'phi'
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
    c_E = '-0.2'
    c_dE_dT = '-0.0005'
    c_Z = '2'
    c_k0 = 1e-15
    boundary = 'solid_liquid_interface'
    phi = 0.0 #'phi'
    temperature = 'T_fluid'
    subdomain1 = ${fluid_blocks}
    subdomain2 = ${solid_blocks}
    wall_cell_is_bulk = true
  []
[]

[FunctorMaterials]
  [diff_coef_species]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'diff_coef_species'
    subdomain_to_prop_value = 'pipe 1e-6
                               flow 1e-2
                               flow-heated 1e-2
                               flow-cooled 1e-2'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -sub_pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'
  automatic_scaling = true
  nl_rel_tol = 1e-1
  nl_abs_tol = 1e-12
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
