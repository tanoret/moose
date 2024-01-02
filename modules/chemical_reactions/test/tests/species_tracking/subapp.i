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
  [something]
    type = MagicBookSpecies
    isotopes_name = 'H1 H2 H3 Be4 Be5'
    # initial_value_isotopes = 2.0
    # intial_values_isotopes = '2.0 2.0 2.0 2.0 2.0'
    uo_function = 'isotopes_to_elements'
    sub_filenames = 'subapp.i'
  []
[]

[Postprocessors]
  [average_h1_sub]
    type = ElementAverageValue
    variable = H1
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
  [csv]
    type = CSV
    file_base = 'sub_app'
    time_data = true
  []
[]
