n1 = 1.46
n2 = 1.00

alpha = 0.719062863767435
alpha1 = 0.2115609606560856
alpha2 = 0.5970451183268032
beta1 = 0.19452369229011088
beta2 = -0.4203259177124563
eta1 = -2.623417821661131
eta2 = 9.947147040979628

# alpha = 0.92
# alpha1 = 2.398432097141336
# alpha2 = 1.1432345695253308
# beta1 = 0.04711223359321187
# beta2 = 0.16122109974012147
# eta1 = 32.16555071624166
# eta2 = 14.958339087605234

k = 1.0
h = 1.0
Tb = 300
T0 = 1000
epsilon = 0.1

# k_rhocp = 5.5e-7

nu_min = 1e-2
nu1 = 4.28e13
nu_high = 1e16

[Mesh]
    [cmg]
        type = CartesianMeshGenerator
        dim = 1
        dx = '1'
        ix = '100'
    []
[]

[Functions]
  [kappa]
    type = ParsedFunction
    expression = "if((x<0.5), 0.1, 10.0)"
  []
[]

[Variables]
  [T]
    type = MooseVariableFVReal
    initial_condition = ${T0}
  []

  [psi1]
    type = MooseVariableFVReal
  []

  [psi2]
    type = MooseVariableFVReal
  []
[]

[FVKernels]
  [diffusion1]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi1
    epsilon = ${epsilon}
    kappa = kappa
    order = first
    # boundaries_to_avoid = 'left right'
    # force_boundary_execution = true
  []

  [diffusion2]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi2
    epsilon = ${epsilon}
    kappa = kappa
    order = second
    # boundaries_to_avoid = 'left right'
    # force_boundary_execution = true
  []

  # [source1]
  #   type = FVSP3ThermalRadiationSourceSink
  #   variable = psi1
  #   T = 'T'
  #   nu = ${nu1}
  #   nu_low =${nu1}
  #   nu_high = ${nu_high}
  #   refraction_index = ${n1}
  #   kappa = kappa
  #   # boundaries_to_avoid = 'left right'
  # []

  # [source2]
  #   type = FVSP3ThermalRadiationSourceSink
  #   variable = psi2
  #   T = 'T'
  #   nu = ${nu1}
  #   nu_low = ${nu1}
  #   nu_high = ${nu_high}
  #   refraction_index = ${n1}
  #   kappa = kappa
  #   # boundaries_to_avoid = 'left right'
  # []

  [sink1]
    type = FVSP3ThermalRadiationSink
    variable = psi1
    T = 'T'
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu_high}
    refraction_index = ${n1}
    kappa = kappa
    # boundaries_to_avoid = 'left right'
  []

  [sink2]
    type = FVSP3ThermalRadiationSink
    variable = psi2
    T = 'T'
    nu = ${nu1}
    nu_low = ${nu1}
    nu_high = ${nu_high}
    refraction_index = ${n1}
    kappa = kappa
    # boundaries_to_avoid = 'left right'
  []

  [source1]
    type = FVSP3ThermalRadiationSource
    variable = psi1
    T = 'T'
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu_high}
    refraction_index = ${n1}
    kappa = kappa
    # boundaries_to_avoid = 'left right'
  []

  [source2]
    type = FVSP3ThermalRadiationSource
    variable = psi2
    T = 'T'
    nu = ${nu1}
    nu_low = ${nu1}
    nu_high = ${nu_high}
    refraction_index = ${n1}
    kappa = kappa
    # boundaries_to_avoid = 'left right'
  []

  [energy_source]
    type = FVSP3TemperatureSourceSink
    variable = T
    absorptivities = 'kappa'
    psi_1 = 'psi1'
    psi_2 = 'psi2'
    force_boundary_execution = true
  []

  [energy_time]
    type = FVTimeKernel
    variable = T
  []

  [energy_diffusion]
    type = FVDiffusion
    variable = T
    coeff = ${k}
  []
[]

[FVBCs]
  [BC1]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right'
    variable = psi1
    Tb = ${Tb}
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu_high}
    refraction_index = ${n1}
    kappa = kappa
    epsilon = ${epsilon}
    psi = 'psi2'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC2]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right'
    variable = psi2
    Tb = ${Tb}
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu_high}
    refraction_index = ${n1}
    kappa = kappa
    epsilon = ${epsilon}
    psi = 'psi1'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC_temp_N]
    type = FVNeumannBC
    boundary = 'left'
    variable = T
    value = 0
  []

  [BC_temperature]
    # type = FVFunctorConvectiveHeatFluxBC
    # boundary = 'right'
    # variable = T
    # T_bulk = ${Tb}
    # T_solid = T
    # is_solid = true
    # heat_transfer_coefficient = ${h}

    type = FVSP3TemperatureBC
    boundary = 'left right'
    variable = T
    Tb = ${Tb}
    n1 = ${n1}
    n2 = ${n2}
    h = ${h}
    k = ${k}
    epsilon = ${epsilon}
    alpha = ${alpha}
    nu1 = ${nu1}
    nu_min = ${nu_min}
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
  dt = 2e-6
  end_time = 0.01

  [TimeIntegrator]
    type = ImplicitEuler
  []

  l_tol = 1e-8
[]

  
[Outputs]
    exodus = true
[]