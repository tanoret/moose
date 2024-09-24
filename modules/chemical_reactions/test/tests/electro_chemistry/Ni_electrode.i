F = 96485.3321
epsilon_0 = 8.854e-12
epsilon_r = 80
T = 800
E_Ni = -0.23
dE_dT_Ni = 0.146e-3
Foeps = ${fparse F/(epsilon_0*epsilon_r)}

[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 1
    dx = '0.1 1.0'
    ix = '4 10'
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
  [Ni_s]
    type = MooseVariableFVReal
    initial_condition = 1.0
    block = '1'
  []
  [Ni]
    type = MooseVariableFVReal
    initial_condition = 0.1
    block = '2'
  []
  # [phi]
  #   type = MooseVariableFVReal
  #   initial_condition = 0.0
  #   block = '1 2'
  # []
[]

[FVKernels]
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

  # [laplace_phi]
  #   type = FVDiffusion
  #   variable = 'phi'
  #   coeff = '1.0'
  #   block = '1 2'
  # []
  # [Ni_source]
  #   type = FVCoupledForce
  #   variable = 'phi'
  #   v = 'Ni_ex'
  #   coef = '${Foeps}'
  #   block = '2'
  # []
[]

[FVInterfaceKernels]
  [interface_Ni]
    type = FVButlerBolmerInterface
    variable1 = 'Ni'
    variable2 = 'Ni_s'
    c = 'Ni'
    c_solid = 'Ni_s'
    c_E = ${E_Ni}
    c_dE_dT = ${dE_dT_Ni}
    c_Z = '2'
    c_k0 = 1.0
    T_ref = 300.0
    boundary = 'solid_liquid_interface'
    phi = '0' #'phi'
    temperature = ${T}
    subdomain1 = '2'
    subdomain2 = '1'
    wall_cell_is_bulk = true
  []
[]

[FVBCs]
  # [phi]
  #   type = FVDirichletBC
  #   variable = 'phi'
  #   boundary = 'left right'
  #   value = 0.0
  # []
[]

[FunctorMaterials]
  [diff_coef]
    type = PiecewiseByBlockFunctorMaterial
    prop_name = 'diff_coef'
    subdomain_to_prop_value = '1 1.0
                               2 1.0'
  []
  # [Ni_ex]
  #   type = ADParsedFunctorMaterial
  #   property_name = 'Ni_ex'
  #   functor_names = 'Ni'
  #   functor_symbols = 'Ni'
  #   expression = 'Ni - 0.0'
  # []
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
  end_time = 1e10
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  automatic_scaling = true
  # off_diagonals_in_auto_scaling = true
  # compute_scaling_once = false
  # residual_and_jacobian_together = true
  nl_abs_tol = 1e-15
  nl_rel_tol = 1e-6
  nl_max_its = 50
  steady_state_detection = false
  steady_state_tolerance = 1e-12
[]

[Postprocessors]
  [average_Ni_s]
    type = ElementAverageValue
    variable = Ni_s
    block = 1
  []
  [average_Ni_l]
    type = ElementAverageValue
    variable = Ni
    block = 2
  []
  [ratio_Ni_l_o_Ni_s]
    type = ParsedPostprocessor
    expression = 'average_Ni_l / average_Ni_s'
    pp_names = 'average_Ni_l average_Ni_s'
  []
[]

[Outputs]
  exodus = true
[]
