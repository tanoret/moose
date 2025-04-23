k = 1.0
h = 1.0
nu1 = 4.28e13
n1 = 1.46
n2 = 1
alpha = 0.92
# cm = 900
# rhom = 2200
Tb = 300
T0 = 1000

eps = 1.0
kap = 1.0

[Mesh]
    [cmg]
        type = CartesianMeshGenerator
        dim = 1
        dx = '1'
        ix = '100'
    []
[]

[Variables]
    [T]
        type = MooseVariableFVReal
        initial_condition = ${T0}
    []
    [psi_1]
        type = MooseVariableFVReal
        initial_condition = 0.0
    []
    [psi_2]
        type = MooseVariableFVReal
        initial_condition = 0.0
    []
[]

[FVKernels]
    [diffusion_1]
        type = FVSP3ThermalRadiationDiffusion
        variable = psi_1
        epsilon = ${eps}
        kappa = ${kap}
        order = first  
    []
    [diffusion_2]
        type = FVSP3ThermalRadiationDiffusion
        variable = psi_2
        epsilon = ${eps}
        kappa = ${kap}
        order = second 
    []
    [source_1]
        type = FVSP3ThermalRadiationSourceSink
        variable = psi_1
        T = 'T'
        nu = ${nu1}
        refraction_index = ${n1}
        kappa = ${kap}
    []
    [source_2]
        type = FVSP3ThermalRadiationSourceSink
        variable = psi_2
        T = 'T'
        nu = ${nu1}
        refraction_index = ${n1}
        kappa = ${kap}
    []

    [energy_source]
        type = FVSP3TemperatureSourceSink
        variable = T
        absorptivities = '${kap}'
        psi_1 = 'psi_1'
        psi_2 = 'psi_2'
        band_frequency_width = '${nu1}'
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
    [BC_1]
        type = FVSP3ThermalRadiationBC
        boundary = 'left right'
        variable = psi_1
        T = 'T'
        nu = ${nu1}
        refraction_index = ${n1}
        kappa = ${kap}
        epsilon = ${eps}
        psi = 'psi_2'
        order = first
    []
    [BC_2]
        type = FVSP3ThermalRadiationBC
        boundary = 'left right'
        variable = psi_2
        T = 'T'
        nu = ${nu1}
        refraction_index = ${n2}
        kappa = ${kap}
        epsilon = ${eps}
        psi = 'psi_1'
        order = second
    []

    [BC_temperature]
        type = FVSP3TemperatureBC
        boundary = 'left right'
        variable = T
        Tb = ${Tb}
        n1 = ${n1}
        n2 = ${n2}
        h = ${h}
        k = ${k}
        alpha = ${alpha}
        nu1 = ${nu1}
        Nint = 1
    []
[]

[Executioner]
    type = Transient
    steady_state_detection = true
    steady_state_tolerance = 1e-12 
    dt = 1e-2
    num_steps = 100
[]
  
[Outputs]
    exodus = true
[]