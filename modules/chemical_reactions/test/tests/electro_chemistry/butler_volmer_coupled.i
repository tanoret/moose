F = 96485.3321
epsilon_0 = 8.854e-12
epsilon_r = 80
Foeps = ${fparse F/(epsilon_0*epsilon_r)}

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.1 1.0'
    dy = '1.0'
    ix = '3 10'
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
  [Zn_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '1'
  []
  [Zn]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '2'
  []
  [Ni_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '1'
  []
  [Ni]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '2'
  []
  [phi]
    type = MooseVariableFVReal
    initial_condition = 0.0
    block = '1 2'
  []
[]

[FVKernels]

  [diffusion_solid_Zn]
    type = FVDiffusion
    variable = 'Zn_s'
    coeff = 'diff_coef'
    block = '1'
  []
  [diffusion_liquid_Zn]
    type = FVDiffusion
    variable = 'Zn'
    coeff = 'diff_coef'
    block = '2'
  []

  [diffusion_solid_Ni]
    type = FVDiffusion
    variable = 'Ni_s'
    coeff = 'diff_coef'
    block = '1'
  []
  [diffusion_liquid_Ni]
    type = FVDiffusion
    variable = 'Ni'
    coeff = 'diff_coef'
    block = '2'
  []

  [laplace_phi]
    type = FVDiffusion
    variable = 'phi'
    coeff = '1.0'
    block = '1 2'
  []
  [Ni_source]
    type = FVCoupledForce
    variable = 'phi'
    v = 'Ni'
    coef = '${Foeps}'
  []
  [Zn_source]
    type = FVCoupledForce
    variable = 'phi'
    v = 'Zn'
    coef = '${fparse 2.0*Foeps}'
  []
[]

[FVInterfaceKernels]
  [interface_Zn]
    type = FVButlerBolmerInterface
    variable1 = 'Zn'
    variable2 = 'Zn_s'
    c = 'Zn'
    c_solid = 'Zn_s'
    c_E = '-0.76'
    c_Z = '2'
    c_k0 = 1e-6
    boundary = 'solid_liquid_interface'
    phi = 'phi'
    temperature = 300.0
    subdomain1 = '2'
    subdomain2 = '1'
    wall_cell_is_bulk = true
  []
  [interface_Ni]
    type = FVButlerBolmerInterface
    variable1 = 'Ni'
    variable2 = 'Ni_s'
    c = 'Ni'
    c_solid = 'Ni_s'
    c_E = '-0.23'
    c_Z = '2'
    c_k0 = 1e-6
    boundary = 'solid_liquid_interface'
    phi = 'phi'
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
    value = 0.0
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
    value = 0.0
  []

  [phi]
    type = FVDirichletBC
    variable = 'phi'
    boundary = 'left right'
    value = 0.0
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
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  nl_abs_tol = 1e-12
  nl_rel_tol = 1e-12
  nl_max_its = 50
  residual_and_jacobian_together = true
[]

[Outputs]
  exodus = true
[]
