sigma_0 = 1.0
sigma_tr_1 = 1.0
source = 1.0

[Mesh]
    [cmg]
        type = CartesianMeshGenerator
        dim = 2
        dx = '1.0'
        dy = '1.0'
        ix = '20'
        iy = '20'
    []
[]

[Variables]
    [phi_0]
        type = MooseVariableFVReal
    []
[]

[FVKernels]
    [diffusion]
        type = FVPNThermalRadiation
        variable = phi_0
        n = 0
        sigma_n_plus_1 = ${sigma_tr_1}
    []
    [reaction]
        type = FVReaction
        variable = phi_0
        rate = ${sigma_0}
    []
    [volumetric_source]
        type = FVBodyForce
        variable = phi_0
        function = ${source}
    []
[]

[FVBCs]
    # Top and bottom are reflective. So, no BCs
    [left_BC]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'left'
        phi_0 = 'phi_0'
        n = 0
    []
    [right_BC]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'right'
        phi_0 = 'phi_0'
        n = 0
    []
[]

[Executioner]
    type = Steady
[]

[Outputs]
    exodus = true
[]