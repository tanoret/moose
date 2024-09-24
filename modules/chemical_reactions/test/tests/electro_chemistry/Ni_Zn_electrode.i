T = 800
T_ref = 300
E_Ni = -0.23
dE_dT_Ni = 0.146e-3
E_Zn = -0.76
dE_dT_Zn = 0.119e-3
phi_back_Ni = ${fparse -(E_Ni + dE_dT_Ni*(T-T_ref))}
phi_back_Zn = ${fparse -(E_Zn + dE_dT_Zn*(T-T_ref))}

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 1
    dx = '0.1 1.0'
    ix = '10 10'
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
  [Zn_time]
    type = FVTimeKernel
    variable = 'Zn'
  []
  [diffusion_liquid_Zn]
    type = FVDiffusion
    variable = 'Zn'
    coeff = 'diff_coef'
    block = '2'
  []
  [Zn_s_time]
    type = FVTimeKernel
    variable = 'Zn_s'
  []
  [diffusion_solid_Zn]
    type = FVDiffusion
    variable = 'Zn_s'
    coeff = 'diff_coef'
    block = '1'
  []
  [Ni_time]
    type = FVTimeKernel
    variable = 'Ni'
  []
  [diffusion_liquid_Ni]
    type = FVDiffusion
    variable = 'Ni'
    coeff = 'diff_coef'
    block = '2'
  []
  [Ni_s_time]
    type = FVTimeKernel
    variable = 'Ni_s'
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
    variable1 = 'Zn'
    variable2 = 'Zn_s'
    c = 'Zn'
    c_solid = 'Zn_s'
    c_E = ${E_Zn}
    c_dE_dT = ${dE_dT_Zn}
    T_ref = ${T_ref}
    c_Z = '2'
    c_k0 = 1.0
    boundary = 'solid_liquid_interface'
    phi = ${phi_back_Ni} #'phi'
    temperature = ${T}
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
    c_E = ${E_Ni}
    c_dE_dT = ${dE_dT_Ni}
    T_ref = ${T_ref}
    c_Z = '2'
    c_k0 = 1.0
    boundary = 'solid_liquid_interface'
    phi = ${phi_back_Zn}  #'phi'
    temperature = ${T}
    subdomain1 = '2'
    subdomain2 = '1'
    wall_cell_is_bulk = true
  []
[]

[FunctorMaterials]
  [diff_coef]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'diff_coef'
    subdomain_to_prop_value = '1 1.0
                               2 1.0'
  []
[]

[Preconditioning]
  [smp]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  start_time = 0.0
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 0.1
    optimal_iterations = 100
    iteration_window = 2
  []
  end_time = 1e12
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  automatic_scaling = true
  nl_abs_tol = 1e-14
  nl_rel_tol = 1e-6
  nl_max_its = 50
  steady_state_detection = false
  steady_state_tolerance = 1e-12
[]

[Postprocessors]
  [average_Zn_l]
    type = ElementAverageValue
    variable = Zn
    block = 2
  []
  [average_Ni_l]
    type = ElementAverageValue
    variable = Ni
    block = 2
  []
  [average_Zn_s]
    type = ElementAverageValue
    variable = Zn_s
    block = 1
  []
  [average_Ni_s]
    type = ElementAverageValue
    variable = Ni_s
    block = 1
  []
  [ratio_Ni_l_o_Ni_s]
    type = ParsedPostprocessor
    expression = 'average_Ni_l / average_Ni_s'
    pp_names = 'average_Ni_l average_Ni_s'
  []
  [ratio_Zn_l_o_Zn_s]
    type = ParsedPostprocessor
    expression = 'average_Zn_l / average_Zn_s'
    pp_names = 'average_Zn_l average_Zn_s'
  []
[]

[Outputs]
  exodus = true
[]
