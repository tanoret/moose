[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 3
    nx = 2
    ny = 2
    nz = 2
  []
[]

[AuxVariables]
  # [H1]
  #   type = MooseVariableFVReal
  # []
  # [H2]
  #   type = MooseVariableFVReal
  # []
  # [Be4]
  #   type = MooseVariableFVReal
  # []
  # [Be5]
  #   type = MooseVariableFVReal
  # []
[]

[UserObjects]
  [isotopes_tracker]
    type = MagicBookSpecies
    isotopes_name = 'H1 H2 H3 Be4 Be5'
    # initial_value_isotopes = 2.0
    intial_values_isotopes = '2.0 2.0 2.0 2.0 2.0'
    uo_function = 'elements_to_isotopes'
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
