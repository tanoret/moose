sigma_0 = 1.0
sigma_2 = 1.0
sigma_tr_1 = 1.0
sigma_tr_3 = 1.0
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
    [phi_2]
        type = MooseVariableFVReal
    []
[]

[FVKernels]
    [diffusion_0]
        type = FVPNThermalRadiation
        variable = phi_0
        n = 0
        sigma_n_plus_1 = ${sigma_tr_1}
        phi_n_plus_2 = phi_2
    []
    [reaction_0]
        type = FVReaction
        variable = phi_0
        rate = ${sigma_0}
    []
    [volumetric_source_0]
        type = FVBodyForce
        variable = phi_0
        function = ${source}
    []
    [diffusion_2]
        type = FVPNThermalRadiation
        variable = phi_2
        n = 2
        sigma_n_minus_1 = ${sigma_tr_1}
        sigma_n_plus_1 = ${sigma_tr_3}
        phi_n_minus_2 = phi_0
    []
    [reaction_2]
        type = FVReaction
        variable = phi_2
        rate = ${sigma_2}
    []
[]

[FVBCs]
    # Top and bottom are reflective. So, no BCs
    [left_BC_0]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'left'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        n = 0
    []
    [right_BC_0]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'right'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        n = 0
    []
    [left_BC_2]
        type = FVPNMarshakRadiativeBC
        variable = phi_2
        boundary = 'left'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        n = 2
    []
    [right_BC_2]
        type = FVPNMarshakRadiativeBC
        variable = phi_2
        boundary = 'right'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        n = 2
    []
[]

[Executioner]
    type = Steady
[]

[Outputs]
    exodus = true
[]