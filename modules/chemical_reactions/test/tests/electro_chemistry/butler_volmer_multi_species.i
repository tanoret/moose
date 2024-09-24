[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.1 1.0'
    dy = '1.0'
    ix = '10 10'
    iy = '10'
    subdomain_id = '1 2'
  []
  [solid_liquid_interface]
    type = SideSetsBetweenSubdomainsGenerator
    input = cmg
    new_boundary = 'solid_liquid_interface'
    primary_block = '1'
    paired_block = '2'
  []
[]

[Variables]
  [Zn]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '2'
  []
  [Zn_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '1'
  []
  [Ni]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '2'
  []
  [Ni_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '1'
  []
[]

[FVKernels]
  [diffusion_liquid_Zn]
    type = FVDiffusion
    variable = 'Zn'
    coeff = 'diff_coef'
    block = '2'
  []
  [diffusion_solid_Zn]
    type = FVDiffusion
    variable = 'Zn_s'
    coeff = 'diff_coef'
    block = '1'
  []
  [diffusion_liquid_Ni]
    type = FVDiffusion
    variable = 'Ni'
    coeff = 'diff_coef'
    block = '2'
  []
  [diffusion_solid_Ni]
    type = FVDiffusion
    variable = 'Ni_s'
    coeff = 'diff_coef'
    block = '1'
  []
[]

[FVInterfaceKernels]
  [interface_Zn]
    type = FVButlerBolmerInterface
    interface_formulation = 'multi-species'
    variable1 = 'Zn'
    variable2 = 'Zn_s'
    c = 'Zn'
    c_solid = 'Zn_s'
    c_E = '-0.76'
    c_Z = '2'
    coupled_c = 'Ni'
    coupled_c_solid = 'Ni_s'
    coupled_E = '-0.23'
    coupled_Z = '2'
    coupled_k0 = '1e-10'
    boundary = 'solid_liquid_interface'
    phi = '0.0'
    temperature = 300.0
    subdomain1 = '2'
    subdomain2 = '1'
    wall_cell_is_bulk = true
  []
  [interface_Ni]
    type = FVButlerBolmerInterface
    interface_formulation = 'multi-species'
    variable1 = 'Ni'
    variable2 = 'Ni_s'
    c = 'Ni'
    c_solid = 'Ni_s'
    c_E = '-0.23'
    c_Z = '2'
    coupled_c = 'Zn'
    coupled_c_solid = 'Zn_s'
    coupled_E = '-0.76'
    coupled_Z = '2'
    coupled_k0 = '1e-10'
    boundary = 'solid_liquid_interface'
    phi = '0.0'
    temperature = 300.0
    subdomain1 = '2'
    subdomain2 = '1'
    wall_cell_is_bulk = true
  []
[]

[FVBCs]
  [left_Zn]
    type = FVDirichletBC
    variable = 'Zn_s'
    boundary = 'left'
    value = 1.0
  []

  [right_Zn]
    type = FVDirichletBC
    variable = 'Zn'
    boundary = 'right'
    value = 1.0
  []

  [left_Ni]
    type = FVDirichletBC
    variable = 'Ni_s'
    boundary = 'left'
    value = 1.0
  []

  [right_Ni]
    type = FVDirichletBC
    variable = 'Ni'
    boundary = 'right'
    value = 1.0
  []
[]

[FunctorMaterials]
  [diff_coef]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'diff_coef'
    subdomain_to_prop_value = '1 1e-10
                               2 1e-3'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  # end_time = 1.0
  # dt = 0.1
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-12
  nl_max_its = 50
  residual_and_jacobian_together = true
[]

[Outputs]
  exodus = true
[]
