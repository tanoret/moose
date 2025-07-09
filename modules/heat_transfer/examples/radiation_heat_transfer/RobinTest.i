k = 2
h = 5
Q = 100
Tb = 300
T0 = 300
alpha = 0.0
epsilon = 1.0

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 1
    xmin = 0.0
    xmax = 2.0
    nx = 100
  []
[]

[Variables]
  [T]
    type = MooseVariableFVReal
    initial_condition = ${T0}
  []
[]

[FVKernels]
  [energy_source]
    type = FVBodyForce
    variable = T
    value = ${Q}
  []

  # [energy_time]
  #   type = FVTimeKernel
  #   variable = T
  # []

  [energy_diffusion]
    type = FVDiffusion
    variable = T
    coeff = ${k}
  []
[]

[FVBCs]
  [BC_temperature]
    type = FVSP3TemperatureBC
    boundary = 'right'
    variable = T
    Tb = ${Tb}
    n1 = 1.0 #dummy
    n2 = 1.0 #dummy
    h = ${h}
    k = ${k}
    epsilon = ${epsilon}
    alpha = ${alpha}
    nu1 = 1.0    #dummy
    nu_min = 0.1 #dummy
  []

  # [rightBC]
  #   type = FVDirichletBC
  #   boundary = 'right'
  #   value = ${Tb}
  #   variable = T
  # []

  [leftBC]
  type = FVNeumannBC
    boundary = 'left' 
    value = 0.0
    variable = T
  []
[]

[Executioner]
  type = Transient
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       NONZERO                superlu_dist' 
  # petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'
  # petsc_options_value = 'lu   superlu_dist' 

  nl_rel_tol = 5e-6
  nl_abs_tol = 5e-6

  start_time = 0.0
  dt = 0.1
  end_time = 0.1

  [TimeIntegrator]
    type = ImplicitEuler
  []

  l_tol = 1e-8
[]

  
[Outputs]
    exodus = true
[]