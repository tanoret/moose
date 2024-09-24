[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 1.
    ymin = 0
    ymax = 1.
    nx = 20
    ny = 20
  []
[]

[Variables]
  [c]
    type = MooseVariableFVReal
  []
[]

[AuxVariables]
  [phi]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [electric_potential]
    type = ParsedAux
    variable = phi
    use_xyzt = true
    expression = 'x'
  []
[]

[FVKernels]
  [electrophoresis]
    type = FVElectrophoresisSource
    variable = c
    temperature = 1000.0
    D = 1.0
    phi = phi
  []
  [diffusion]
    type = FVDiffusion
    variable = c
    coeff = 1.0
  []
[]

[FVBCs]
  [left]
    type = FVDirichletBC
    variable = c
    boundary = 'left bottom'
    value = 1
  []

  [right]
    type = FVDirichletBC
    variable = c
    boundary = 'right top'
    value = 0
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
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu NONZERO'
  nl_rel_tol = 1e-12
  residual_and_jacobian_together = true
[]

[Outputs]
  exodus = true
[]
