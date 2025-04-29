k = 1.0
h = 1.0
n1 = 1.46
n2 = 1
alpha = 0.92
Tb = 300
T0 = 1000
epsilon = 1.0

nu1 = 42827494000000.0

nu1half = 46121916615384.62
nu2half = 52137818782608.695
nu3half = 59958491600000.01
nu4half = 70539401882352.94
nu5half = 79944655466666.67
nu6half = 92243833230769.23
nu7half = 187370286250000.0
nu8half = 85654988000000.0

kappa1 = 7136.06
kappa2 = 567.32
kappa3 = 267.98
kappa4 = 27.98
kappa5 = 15.45
kappa6 = 7.7
kappa7 = 0.5
kappa8 = 0.4

nu1to2 = 7137915666666.664
nu2to3 = 4542309969696.977
nu3to4 = 12112826585858.578
nu4to5 = 8327568277777.781
nu5to6 = 10706873500000.0
nu6to7 = 14275831333333.328
nu7to8 = 1399031470666666.8
nu8to9 = 13990314706666668

[Mesh]
    [ccmg]
      type = ConcentricCircleMeshGenerator
      has_outer_square = false
      radii = '0.41 0.46  0.5'
      num_sectors = 20
      rings = '1 1 1'
      preserve_volumes = false
      smoothing_max_it = 2
    []
[]

[Variables]
  [T]
    type = MooseVariableFVReal
    initial_condition = ${T0}
  []

  [psi11]
      type = MooseVariableFVReal
  []

  [psi21]
      type = MooseVariableFVReal
  []

  [psi12]
      type = MooseVariableFVReal
  []

  [psi22]
      type = MooseVariableFVReal
  []

  [psi13]
      type = MooseVariableFVReal
  []

  [psi23]
      type = MooseVariableFVReal
  []

  [psi14]
      type = MooseVariableFVReal
  []

  [psi24]
      type = MooseVariableFVReal
  []

  [psi15]
      type = MooseVariableFVReal
  []

  [psi25]
      type = MooseVariableFVReal
  []

  [psi16]
      type = MooseVariableFVReal
  []

  [psi26]
      type = MooseVariableFVReal
  []

  [psi17]
      type = MooseVariableFVReal
  []

  [psi27]
      type = MooseVariableFVReal
  []

  [psi18]
      type = MooseVariableFVReal
  []

  [psi28]
      type = MooseVariableFVReal
  []

[]

[FVKernels]
  [diffusion11]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi11
      epsilon = ${epsilon}
      kappa = ${kappa1}
      order = first
  []

  [diffusion21]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi21
      epsilon = ${epsilon}
      kappa = ${kappa1}
      order = second
  []

  [diffusion12]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi12
      epsilon = ${epsilon}
      kappa = ${kappa2}
      order = first
  []

  [diffusion22]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi22
      epsilon = ${epsilon}
      kappa = ${kappa2}
      order = second
  []

  [diffusion13]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi13
      epsilon = ${epsilon}
      kappa = ${kappa3}
      order = first
  []

  [diffusion23]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi23
      epsilon = ${epsilon}
      kappa = ${kappa3}
      order = second
  []

  [diffusion14]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi14
      epsilon = ${epsilon}
      kappa = ${kappa4}
      order = first
  []

  [diffusion24]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi24
      epsilon = ${epsilon}
      kappa = ${kappa4}
      order = second
  []

  [diffusion15]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi15
      epsilon = ${epsilon}
      kappa = ${kappa5}
      order = first
  []

  [diffusion25]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi25
      epsilon = ${epsilon}
      kappa = ${kappa5}
      order = second
  []

  [diffusion16]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi16
      epsilon = ${epsilon}
      kappa = ${kappa6}
      order = first
  []

  [diffusion26]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi26
      epsilon = ${epsilon}
      kappa = ${kappa6}
      order = second
  []

  [diffusion17]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi17
      epsilon = ${epsilon}
      kappa = ${kappa7}
      order = first
  []

  [diffusion27]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi27
      epsilon = ${epsilon}
      kappa = ${kappa7}
      order = second
  []

  [diffusion18]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi18
      epsilon = ${epsilon}
      kappa = ${kappa8}
      order = first
  []

  [diffusion28]
      type = FVSP3ThermalRadiationDiffusion
      variable = psi28
      epsilon = ${epsilon}
      kappa = ${kappa8}
      order = second
  []

  [source11]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi11
      T = 'T'
      nu = ${nu1half}
      refraction_index = ${n1}
      kappa = ${kappa1}
  []

  [source21]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi21
      T = 'T'
      nu = ${nu1half}
      refraction_index = ${n1}
      kappa = ${kappa1}
  []

  [source12]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi12
      T = 'T'
      nu = ${nu2half}
      refraction_index = ${n1}
      kappa = ${kappa2}
  []

  [source22]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi22
      T = 'T'
      nu = ${nu2half}
      refraction_index = ${n1}
      kappa = ${kappa2}
  []

  [source13]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi13
      T = 'T'
      nu = ${nu3half}
      refraction_index = ${n1}
      kappa = ${kappa3}
  []

  [source23]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi23
      T = 'T'
      nu = ${nu3half}
      refraction_index = ${n1}
      kappa = ${kappa3}
  []

  [source14]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi14
      T = 'T'
      nu = ${nu4half}
      refraction_index = ${n1}
      kappa = ${kappa4}
  []

  [source24]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi24
      T = 'T'
      nu = ${nu4half}
      refraction_index = ${n1}
      kappa = ${kappa4}
  []

  [source15]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi15
      T = 'T'
      nu = ${nu5half}
      refraction_index = ${n1}
      kappa = ${kappa5}
  []

  [source25]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi25
      T = 'T'
      nu = ${nu5half}
      refraction_index = ${n1}
      kappa = ${kappa5}
  []

  [source16]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi16
      T = 'T'
      nu = ${nu6half}
      refraction_index = ${n1}
      kappa = ${kappa6}
  []

  [source26]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi26
      T = 'T'
      nu = ${nu6half}
      refraction_index = ${n1}
      kappa = ${kappa6}
  []

  [source17]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi17
      T = 'T'
      nu = ${nu7half}
      refraction_index = ${n1}
      kappa = ${kappa7}
  []

  [source27]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi27
      T = 'T'
      nu = ${nu7half}
      refraction_index = ${n1}
      kappa = ${kappa7}
  []

  [source18]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi18
      T = 'T'
      nu = ${nu8half}
      refraction_index = ${n1}
      kappa = ${kappa8}
  []

  [source28]
      type = FVSP3ThermalRadiationSourceSink
      variable = psi28
      T = 'T'
      nu = ${nu8half}
      refraction_index = ${n1}
      kappa = ${kappa8}
  []

  [energy_source]
      type = FVSP3TemperatureSourceSink
      variable = T
      absorptivities = '${kappa1} ${kappa2} ${kappa3} ${kappa4} ${kappa5} ${kappa6} ${kappa7} ${kappa8} '
      psi_1 = 'psi11 psi12 psi13 psi14 psi15 psi16 psi17 psi18 '
      psi_2 = 'psi21 psi22 psi23 psi24 psi25 psi26 psi27 psi28 '
      band_frequency_width = '${nu1to2} ${nu2to3} ${nu3to4} ${nu4to5} ${nu5to6} ${nu6to7} ${nu7to8} ${nu8to9} '
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
  [BC11]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi11
      T = 'T'
      nu = ${nu1half}
      refraction_index = ${n1}
      kappa = ${kappa1}
      epsilon = ${epsilon}
      psi = 'psi21'
      order = first
  []

  [BC21]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi21
      T = 'T'
      nu = ${nu1half}
      refraction_index = ${n1}
      kappa = ${kappa1}
      epsilon = ${epsilon}
      psi = 'psi11'
      order = second
  []

  [BC12]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi12
      T = 'T'
      nu = ${nu2half}
      refraction_index = ${n1}
      kappa = ${kappa2}
      epsilon = ${epsilon}
      psi = 'psi22'
      order = first
  []

  [BC22]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi22
      T = 'T'
      nu = ${nu2half}
      refraction_index = ${n1}
      kappa = ${kappa2}
      epsilon = ${epsilon}
      psi = 'psi12'
      order = second
  []

  [BC13]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi13
      T = 'T'
      nu = ${nu3half}
      refraction_index = ${n1}
      kappa = ${kappa3}
      epsilon = ${epsilon}
      psi = 'psi23'
      order = first
  []

  [BC23]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi23
      T = 'T'
      nu = ${nu3half}
      refraction_index = ${n1}
      kappa = ${kappa3}
      epsilon = ${epsilon}
      psi = 'psi13'
      order = second
  []

  [BC14]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi14
      T = 'T'
      nu = ${nu4half}
      refraction_index = ${n1}
      kappa = ${kappa4}
      epsilon = ${epsilon}
      psi = 'psi24'
      order = first
  []

  [BC24]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi24
      T = 'T'
      nu = ${nu4half}
      refraction_index = ${n1}
      kappa = ${kappa4}
      epsilon = ${epsilon}
      psi = 'psi14'
      order = second
  []

  [BC15]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi15
      T = 'T'
      nu = ${nu5half}
      refraction_index = ${n1}
      kappa = ${kappa5}
      epsilon = ${epsilon}
      psi = 'psi25'
      order = first
  []

  [BC25]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi25
      T = 'T'
      nu = ${nu5half}
      refraction_index = ${n1}
      kappa = ${kappa5}
      epsilon = ${epsilon}
      psi = 'psi15'
      order = second
  []

  [BC16]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi16
      T = 'T'
      nu = ${nu6half}
      refraction_index = ${n1}
      kappa = ${kappa6}
      epsilon = ${epsilon}
      psi = 'psi26'
      order = first
  []

  [BC26]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi26
      T = 'T'
      nu = ${nu6half}
      refraction_index = ${n1}
      kappa = ${kappa6}
      epsilon = ${epsilon}
      psi = 'psi16'
      order = second
  []

  [BC17]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi17
      T = 'T'
      nu = ${nu7half}
      refraction_index = ${n1}
      kappa = ${kappa7}
      epsilon = ${epsilon}
      psi = 'psi27'
      order = first
  []

  [BC27]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi27
      T = 'T'
      nu = ${nu7half}
      refraction_index = ${n1}
      kappa = ${kappa7}
      epsilon = ${epsilon}
      psi = 'psi17'
      order = second
  []

  [BC18]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi18
      T = 'T'
      nu = ${nu8half}
      refraction_index = ${n1}
      kappa = ${kappa8}
      epsilon = ${epsilon}
      psi = 'psi28'
      order = first
  []

  [BC28]
      type = FVSP3ThermalRadiationBC
      boundary = outer
      variable = psi28
      T = 'T'
      nu = ${nu8half}
      refraction_index = ${n1}
      kappa = ${kappa8}
      epsilon = ${epsilon}
      psi = 'psi18'
      order = second
  []

    [BC_temperature]
        type = FVSP3TemperatureBC
        boundary =  outer
        variable = T
        Tb = ${Tb}
        n1 = ${n1}
        n2 = ${n2}
        h = ${h}
        k = ${k}
        alpha = ${alpha}
        nu1 = ${nu1}
        Nintegral = 100
    []
[]

[Executioner]
  type = Transient
  solve_type = 'Newton'

  num_steps = 100
  dt = 0.00001

  nl_rel_tol = 1e-6
  l_tol = 1e-6

  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
  petsc_options_value = 'lu       NONZERO                superlu_dist'
[]
  
[Outputs]
    exodus = true
[]