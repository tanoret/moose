n1 = 1.46
n2 = 1.00

alpha = 0.914
alpha1 =  0.084075389781548382
alpha2 =  0.872433935441263020
beta1 = -0.08086509944859642
beta2 = -0.29284032175563224
eta1 =  -2.6234175064678195
eta2 = 9.9471473598607409

k = 1
h = 0.001
Tb = 300
T0 = 1000
epsilon = 1.0
endt = 0.001

# nu_min = 3e-3
# nu1 = 2.9334638e13
# nu2 = 3.4223744e13
# nu3 = 3.7334994e13
# nu4 = 4.5631659e13
# nu5 = 5.1335616e13
# nu6 = 5.8669276e13
# nu7 = 6.8447489e13
# nu8 = 102.6712329e13
# nu_max = 1e17

nu_min = 1e-2
nu1 = 2.93e13
nu2 = 3.42e13
nu3 = 3.73e13
nu4 = 4.56e13
nu5 = 5.13e13
nu6 = 5.87e13
nu7 = 6.84e13
nu8 = 102.67e13
nu_max = 1e16

kappa1to2 = 7136.00
kappa2to3 = 576.32
kappa3to4 = 276.98
kappa4to5 = 27.98
kappa5to6 = 15.45
kappa6to7 = 7.7
kappa7to8 = 0.5
kappa8tomax = 0.4

bias_x = 0.9
bias_y = 0.9
bias_z = 0.9

nx = 31
ny = 31
nz = 31

[Mesh]
  # [gmg]
  #   type = GeneratedMeshGenerator
  #   dim = 3
  #   nx = 61
  #   ny = 61
  #   nz = 61
  #   xmin = 0
  #   xmax = 1
  #   ymin = 0
  #   ymax = 1
  #   zmin = 0
  #   zmax = 1
  #   gauss_lobatto_grid = true
  # []

  # [ref]
  #   type = RefineSidesetGenerator
  #   boundaries = 'left right top bottom front back'
  #   refinement = '1    1    1    1    1    1'
  #   boundary_side = 'both both both both both both'
  #   input = gmg
  # []

  # [tri]
  #   type = ElementsToSimplicesConverter
  #   input = ref
  # []

  # [cmg]
  #   type = CartesianMeshGenerator
  #   dim = 3

  #   dx = '0.00030 0.00035 0.00087 0.00356 0.01492 0.04 0.07 0.07 0.6 0.07 0.07 0.04 0.01492 0.00356 0.00087 0.00035 0.00030'
  #   dy = '0.00030 0.00035 0.00087 0.00356 0.01492 0.04 0.07 0.07 0.6 0.07 0.07 0.04 0.01492 0.00356 0.00087 0.00035 0.00030'
  #   dz = '0.00030 0.00035 0.00087 0.00356 0.01492 0.04 0.07 0.07 0.6 0.07 0.07 0.04 0.01492 0.00356 0.00087 0.00035 0.00030'

  #   # dx = '0.04 0.02 0.03 0.04 0.06 0.11 0.4 0.11 0.06 0.04 0.03 0.02 0.04'
  #   # dy = '0.04 0.02 0.03 0.04 0.06 0.11 0.4 0.11 0.06 0.04 0.03 0.02 0.04'
  #   # dz = '0.04 0.02 0.03 0.04 0.06 0.11 0.4 0.11 0.06 0.04 0.03 0.02 0.04'

  #   ix = '4 4 3 4 4 3 4 2 5 2 4 3 4 4 3 4 4'
  #   iy = '4 4 4 4 4 3 4 2 5 2 4 3 4 4 4 4 4'
  #   iz = '4 4 4 4 4 3 4 2 5 2 4 3 4 4 4 4 4'

  #   # ix = '4 2 2 2 2 2 3 2 2 2 2 2 4'
  #   # iy = '4 2 2 2 2 2 3 2 2 2 2 2 4'
  #   # iz = '4 2 2 2 2 2 3 2 2 2 2 2 4'

  #   # ix = '2 1 1 1 1 1 1 1 1 1 1 1 2'
  #   # iy = '2 1 1 1 1 1 1 1 1 1 1 1 2'
  #   # iz = '2 1 1 1 1 1 1 1 1 1 1 1 2'
  # []
  
  [front_top_rignt]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = 0.0
    xmax = 0.5
    ymin = 0.0
    ymax = 0.5
    zmin = 0.0
    zmax = 0.5
    bias_x = ${fparse bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse bias_z}
  []

  [front_top_left]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = -0.5
    xmax = 0.0
    ymin = 0.0
    ymax = 0.5
    zmin = 0.0
    zmax = 0.5
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse bias_z}
  []

  [front_bottom_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = 0.0
    xmax = 0.5
    ymin = -0.5
    ymax = 0.0
    zmin = 0.0
    zmax = 0.5
    bias_x = ${fparse bias_x}
    bias_y = ${fparse 1/bias_y}
    bias_z = ${fparse bias_z}
  []
  [front_bottom_left]
      type = GeneratedMeshGenerator
      dim = 3
      nx = ${nx}
      ny = ${ny}
      nz = ${nz}
      xmin = -0.5
      xmax = 0.0
      ymin = -0.5
      ymax = 0.0
      zmin = 0.0
      zmax = 0.5
      bias_x = ${fparse 1/bias_x}
      bias_y = ${fparse 1/bias_y}
      bias_z = ${fparse bias_z}
  []
  [back_top_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = 0.0
    xmax = 0.5
    ymin = 0.0
    ymax = 0.5
    zmin = -0.5
    zmax = 0.0
    bias_x = ${fparse bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse 1/bias_z}
  []
  [back_top_left]
      type = GeneratedMeshGenerator
      dim = 3
      nx = ${nx}
      ny = ${ny}
      nz = ${nz}
      xmin = -0.5
      xmax = 0.0
      ymin = 0.0
      ymax = 0.5
      zmin = -0.5
      zmax = 0.0
      bias_x = ${fparse 1/bias_x}
      bias_y = ${fparse bias_y}
      bias_z = ${fparse 1/bias_z}
  []
  [back_bottom_right]
      type = GeneratedMeshGenerator
      dim = 3
      nx = ${nx}
      ny = ${ny}
      nz = ${nz}
      xmin = 0.0
      xmax = 0.5
      ymin = -0.5
      ymax = 0.0
      zmin = -0.5
      zmax = 0.0
      bias_x = ${fparse bias_x}
      bias_y = ${fparse 1/bias_y}
      bias_z = ${fparse 1/bias_z}
  []
  [back_bottom_left]
      type = GeneratedMeshGenerator
      dim = 3
      nx = ${nx}
      ny = ${ny}
      nz = ${nz}
      xmin = -0.5
      xmax = 0.0
      ymin = -0.5
      ymax = 0.0
      zmin = -0.5
      zmax = 0.0
      bias_x = ${fparse 1/bias_x}
      bias_y = ${fparse 1/bias_y}
      bias_z = ${fparse 1/bias_z}
  []
  [join_front_top]
      type = StitchedMeshGenerator
      inputs = 'front_top_rignt front_top_left'
      stitch_boundaries_pairs = 'left right'
  []
  [join_front_bottom]
    type = StitchedMeshGenerator
    inputs = 'front_bottom_right front_bottom_left'
    stitch_boundaries_pairs = 'left right'
  []
  [join_back_top]
    type = StitchedMeshGenerator
    inputs = 'back_top_right back_top_left'
    stitch_boundaries_pairs = 'left right'
  []
  [join_back_bottom]
    type = StitchedMeshGenerator
    inputs = 'back_bottom_right back_bottom_left'
    stitch_boundaries_pairs = 'left right'
  []
  [join_front]
    type = StitchedMeshGenerator
    inputs = 'join_front_top join_front_bottom'
    stitch_boundaries_pairs = 'bottom top'
  []
  [join_back]
    type = StitchedMeshGenerator
    inputs = 'join_back_top join_back_bottom'
    stitch_boundaries_pairs = 'bottom top'
  []
  [joint_all]
    type = StitchedMeshGenerator
    inputs = 'join_front join_back'
    stitch_boundaries_pairs = 'back front'
  []
[]

[Problem]
  nl_sys_names = 'T
                  psi11 psi12 psi13 psi14 psi15 psi16 psi17 psi18
                  psi21 psi22 psi23 psi24 psi25 psi26 psi27 psi28'
[]

[Variables]
  [T]
    type = MooseVariableFVReal
    initial_condition = ${T0}
    solver_sys = 'T'
  []

  [psi11]
    type = MooseVariableFVReal
    solver_sys = 'psi11'
  []

  [psi21]
    type = MooseVariableFVReal
    solver_sys = 'psi21'
  []

  [psi12]
    type = MooseVariableFVReal
    solver_sys = 'psi12'
  []

  [psi22]
    type = MooseVariableFVReal
    solver_sys = 'psi22'
  []

  [psi13]
    type = MooseVariableFVReal
    solver_sys = 'psi13'
  []

  [psi23]
    type = MooseVariableFVReal
    solver_sys = 'psi23'
  []

  [psi14]
    type = MooseVariableFVReal
    solver_sys = 'psi14'
  []

  [psi24]
    type = MooseVariableFVReal
    solver_sys = 'psi24'
  []

  [psi15]
    type = MooseVariableFVReal
    solver_sys = 'psi15'
  []

  [psi25]
    type = MooseVariableFVReal
    solver_sys = 'psi25'
  []

  [psi16]
    type = MooseVariableFVReal
    solver_sys = 'psi16'
  []

  [psi26]
    type = MooseVariableFVReal
    solver_sys = 'psi26'
  []

  [psi17]
    type = MooseVariableFVReal
    solver_sys = 'psi17'
  []

  [psi27]
    type = MooseVariableFVReal
    solver_sys = 'psi27'
  []

  [psi18]
    type = MooseVariableFVReal
    solver_sys = 'psi18'
  []

  [psi28]
    type = MooseVariableFVReal
    solver_sys = 'psi28'
  []
[]

[FVKernels]
  [diffusion11]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi11
    epsilon = ${epsilon}
    kappa = ${kappa1to2}
    order = first
  []

  [diffusion21]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi21
    epsilon = ${epsilon}
    kappa = ${kappa1to2}
    order = second
  []

  [diffusion12]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi12
    epsilon = ${epsilon}
    kappa = ${kappa2to3}
    order = first
  []

  [diffusion22]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi22
    epsilon = ${epsilon}
    kappa = ${kappa2to3}
    order = second
  []

  [diffusion13]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi13
    epsilon = ${epsilon}
    kappa = ${kappa3to4}
    order = first
  []

  [diffusion23]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi23
    epsilon = ${epsilon}
    kappa = ${kappa3to4}
    order = second
  []

  [diffusion14]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi14
    epsilon = ${epsilon}
    kappa = ${kappa4to5}
    order = first
  []

  [diffusion24]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi24
    epsilon = ${epsilon}
    kappa = ${kappa4to5}
    order = second
  []

  [diffusion15]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi15
    epsilon = ${epsilon}
    kappa = ${kappa5to6}
    order = first
  []

  [diffusion25]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi25
    epsilon = ${epsilon}
    kappa = ${kappa5to6}
    order = second
  []

  [diffusion16]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi16
    epsilon = ${epsilon}
    kappa = ${kappa6to7}
    order = first
  []

  [diffusion26]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi26
    epsilon = ${epsilon}
    kappa = ${kappa6to7}
    order = second
  []

  [diffusion17]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi17
    epsilon = ${epsilon}
    kappa = ${kappa7to8}
    order = first
  []

  [diffusion27]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi27
    epsilon = ${epsilon}
    kappa = ${kappa7to8}
    order = second
  []

  [diffusion18]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi18
    epsilon = ${epsilon}
    kappa = ${kappa8tomax}
    order = first
  []

  [diffusion28]
    type = FVSP3ThermalRadiationDiffusion
    variable = psi28
    epsilon = ${epsilon}
    kappa = ${kappa8tomax}
    order = second
  []

  [sink11]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi11
    T = 'T'
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu2}
    refraction_index = ${n1}
    kappa = ${kappa1to2}
  []

  [sink21]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi21
    T = 'T'
    nu = ${nu1}
    nu_low = ${nu1}
    nu_high = ${nu2}
    refraction_index = ${n1}
    kappa = ${kappa1to2}
  []

  [sink12]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi12
    T = 'T'
    nu = ${nu2}
    nu_low =${nu2}
    nu_high = ${nu3}
    refraction_index = ${n1}
    kappa = ${kappa2to3}
  []

  [sink22]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi22
    T = 'T'
    nu = ${nu2}
    nu_low = ${nu2}
    nu_high = ${nu3}
    refraction_index = ${n1}
    kappa = ${kappa2to3}
  []

  [sink13]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi13
    T = 'T'
    nu = ${nu3}
    nu_low =${nu3}
    nu_high = ${nu4}
    refraction_index = ${n1}
    kappa = ${kappa3to4}
  []

  [sink23]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi23
    T = 'T'
    nu = ${nu3}
    nu_low = ${nu3}
    nu_high = ${nu4}
    refraction_index = ${n1}
    kappa = ${kappa3to4}
  []

  [sink14]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi14
    T = 'T'
    nu = ${nu4}
    nu_low =${nu4}
    nu_high = ${nu5}
    refraction_index = ${n1}
    kappa = ${kappa4to5}
  []

  [sink24]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi24
    T = 'T'
    nu = ${nu4}
    nu_low = ${nu4}
    nu_high = ${nu5}
    refraction_index = ${n1}
    kappa = ${kappa4to5}
  []

  [sink15]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi15
    T = 'T'
    nu = ${nu5}
    nu_low =${nu5}
    nu_high = ${nu6}
    refraction_index = ${n1}
    kappa = ${kappa5to6}
  []

  [sink25]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi25
    T = 'T'
    nu = ${nu5}
    nu_low = ${nu5}
    nu_high = ${nu6}
    refraction_index = ${n1}
    kappa = ${kappa5to6}
  []

  [sink16]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi16
    T = 'T'
    nu = ${nu6}
    nu_low =${nu6}
    nu_high = ${nu7}
    refraction_index = ${n1}
    kappa = ${kappa6to7}
  []

  [sink26]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi26
    T = 'T'
    nu = ${nu6}
    nu_low = ${nu6}
    nu_high = ${nu7}
    refraction_index = ${n1}
    kappa = ${kappa6to7}
  []

  [sink17]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi17
    T = 'T'
    nu = ${nu7}
    nu_low =${nu7}
    nu_high = ${nu8}
    refraction_index = ${n1}
    kappa = ${kappa7to8}
  []

  [sink27]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi27
    T = 'T'
    nu = ${nu7}
    nu_low = ${nu7}
    nu_high = ${nu8}
    refraction_index = ${n1}
    kappa = ${kappa7to8}
  []

  [sink18]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi18
    T = 'T'
    nu = ${nu8}
    nu_low =${nu8}
    nu_high = ${nu_max}
    refraction_index = ${n1}
    kappa = ${kappa8tomax}
  []

  [sink28]
    type = FVSP3ThermalRadiationSourceSink
    variable = psi23
    T = 'T'
    nu = ${nu8}
    nu_low = ${nu8}
    nu_high = ${nu_max}
    refraction_index = ${n1}
    kappa = ${kappa8tomax}
  []

  [energy_source]
    type = FVSP3TemperatureSourceSink
    variable = T
    absorptivities = '${kappa1to2} ${kappa2to3} ${kappa3to4} ${kappa4to5} ${kappa5to6} ${kappa6to7} ${kappa7to8} ${kappa8tomax}'
    psi_1 = 'psi11 psi12 psi13 psi14 psi15 psi16 psi17 psi18'
    psi_2 = 'psi21 psi22 psi23 psi24 psi25 psi26 psi27 psi28'
    # force_boundary_execution = true
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
    boundary = 'left right top bottom front back'
    variable = psi11
    Tb = ${Tb}
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu2}
    refraction_index = ${n1}
    kappa = ${kappa1to2}
    epsilon = ${epsilon}
    psi = 'psi21'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC21]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi21
    Tb = ${Tb}
    nu = ${nu1}
    nu_low =${nu1}
    nu_high = ${nu2}
    refraction_index = ${n1}
    kappa = ${kappa1to2}
    epsilon = ${epsilon}
    psi = 'psi11'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC12]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi12
    Tb = ${Tb}
    nu = ${nu2}
    nu_low =${nu2}
    nu_high = ${nu3}
    refraction_index = ${n1}
    kappa = ${kappa2to3}
    epsilon = ${epsilon}
    psi = 'psi22'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC22]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi22
    Tb = ${Tb}
    nu = ${nu2}
    nu_low =${nu2}
    nu_high = ${nu3}
    refraction_index = ${n1}
    kappa = ${kappa2to3}
    epsilon = ${epsilon}
    psi = 'psi12'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC13]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi13
    Tb = ${Tb}
    nu = ${nu3}
    nu_low =${nu3}
    nu_high = ${nu4}
    refraction_index = ${n1}
    kappa = ${kappa3to4}
    epsilon = ${epsilon}
    psi = 'psi23'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC23]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi23
    Tb = ${Tb}
    nu = ${nu3}
    nu_low =${nu3}
    nu_high = ${nu4}
    refraction_index = ${n1}
    kappa = ${kappa3to4}
    epsilon = ${epsilon}
    psi = 'psi13'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC14]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi14
    Tb = ${Tb}
    nu = ${nu4}
    nu_low =${nu4}
    nu_high = ${nu5}
    refraction_index = ${n1}
    kappa = ${kappa4to5}
    epsilon = ${epsilon}
    psi = 'psi24'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC24]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi24
    Tb = ${Tb}
    nu = ${nu4}
    nu_low =${nu4}
    nu_high = ${nu5}
    refraction_index = ${n1}
    kappa = ${kappa4to5}
    epsilon = ${epsilon}
    psi = 'psi14'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC15]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi15
    Tb = ${Tb}
    nu = ${nu5}
    nu_low =${nu5}
    nu_high = ${nu6}
    refraction_index = ${n1}
    kappa = ${kappa5to6}
    epsilon = ${epsilon}
    psi = 'psi25'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC25]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi25
    Tb = ${Tb}
    nu = ${nu5}
    nu_low =${nu5}
    nu_high = ${nu6}
    refraction_index = ${n1}
    kappa = ${kappa5to6}
    epsilon = ${epsilon}
    psi = 'psi15'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC16]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi16
    Tb = ${Tb}
    nu = ${nu6}
    nu_low =${nu6}
    nu_high = ${nu7}
    refraction_index = ${n1}
    kappa = ${kappa6to7}
    epsilon = ${epsilon}
    psi = 'psi26'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC26]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi26
    Tb = ${Tb}
    nu = ${nu6}
    nu_low =${nu6}
    nu_high = ${nu7}
    refraction_index = ${n1}
    kappa = ${kappa6to7}
    epsilon = ${epsilon}
    psi = 'psi16'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC17]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi17
    Tb = ${Tb}
    nu = ${nu7}
    nu_low =${nu7}
    nu_high = ${nu8}
    refraction_index = ${n1}
    kappa = ${kappa7to8}
    epsilon = ${epsilon}
    psi = 'psi27'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC27]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi27
    Tb = ${Tb}
    nu = ${nu7}
    nu_low =${nu7}
    nu_high = ${nu8}
    refraction_index = ${n1}
    kappa = ${kappa7to8}
    epsilon = ${epsilon}
    psi = 'psi17'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC18]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi18
    Tb = ${Tb}
    nu = ${nu8}
    nu_low =${nu8}
    nu_high = ${nu_max}
    refraction_index = ${n1}
    kappa = ${kappa8tomax}
    epsilon = ${epsilon}
    psi = 'psi28'
    order = first
    alpha = ${alpha1}
    beta = ${beta2}
    eta = ${eta1}
  []

  [BC28]
    type = FVSP3ThermalRadiationBC
    boundary = 'left right top bottom front back'
    variable = psi28
    Tb = ${Tb}
    nu = ${nu8}
    nu_low =${nu8}
    nu_high = ${nu_max}
    refraction_index = ${n1}
    kappa = ${kappa8tomax}
    epsilon = ${epsilon}
    psi = 'psi18'
    order = second
    alpha = ${alpha2}
    beta = ${beta1}
    eta = ${eta2}
  []

  [BC_temperature]
    type = FVSP3TemperatureBC
    boundary = 'left right top bottom front back'
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
  # petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_mat_solver_type'
  # petsc_options_value = 'lu       NONZERO                superlu_dist' 
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10000'

  nl_rel_tol = 5e-3
  l_tol = 1e-2
  l_max_its = 100
  nl_max_its = 10000
  
  start_time = 0.0
  dt = 0.00005 #2.5e-6
  end_time = ${endt}
[]

# [VectorPostprocessors]
#   [./lv1]
#     num_points = 60
#     start_point = '0.0 0.5 0.5'
#     end_point = '1.0 0.5 0.5'
#     sort_by = 'x'
#     variable = T
#     type = LineValueSampler
#   [../]
# []

[Outputs]
  exodus = true
  [csv]
    type = CSV
    file_base = Larsen_adaptive3D_
    execute_on = final
  []
[]