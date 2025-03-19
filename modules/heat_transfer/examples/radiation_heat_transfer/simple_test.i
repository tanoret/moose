sigma_0 = 0.5
sigma_2 = 1.0
sigma_4 = 1.0
sigma_tr_1 = 1.0
sigma_tr_3 = 1.0
sigma_tr_5 = 1.0
source = 12.566

[Mesh]
    [cmg]
        type = CartesianMeshGenerator
        dim = 2
        dx = '10'
        dy = '10'
        ix = '40'
        iy = '40'
    []
[]

[Variables]
    [phi_0]
        type = MooseVariableFVReal
    []
    [phi_2]
        type = MooseVariableFVReal
    []
    [phi_4]
        type = MooseVariableFVReal
    []
[]

[Functions]
    [source_fn] 
        type = ParsedFunction
        expression = "if(x>=2 & x<=8 & y>=2 & y<=8, ${source}, 0)"
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
        function = source_fn
    []
    [diffusion_2]
        type = FVPNThermalRadiation
        variable = phi_2
        n = 2
        sigma_n_minus_1 = ${sigma_tr_1}
        sigma_n_plus_1 = ${sigma_tr_3}
        phi_n_minus_2 = phi_0
        phi_n_plus_2 = phi_4
    []
    [reaction_2]
        type = FVReaction
        variable = phi_2
        rate = ${sigma_2}
    []
    [diffusion_4]
        type = FVPNThermalRadiation
        variable = phi_4
        n = 4
        sigma_n_minus_1 = ${sigma_tr_3}
        sigma_n_plus_1 = ${sigma_tr_5}
        phi_n_minus_2 = phi_2
    []
    [reaction_4]
        type = FVReaction
        variable = phi_4
        rate = ${sigma_4}
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
    [top_BC_0]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'top'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        n = 0
    []
    [bottom_BC_0]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'bottom'
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
        phi_4 = 'phi_4'
        n = 2
    []
    [right_BC_2]
        type = FVPNMarshakRadiativeBC
        variable = phi_2
        boundary = 'right'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 2
    []
    [top_BC_2]
        type = FVPNMarshakRadiativeBC
        variable = phi_2
        boundary = 'top'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 2
    []
    [bottom_BC_2]
        type = FVPNMarshakRadiativeBC
        variable = phi_2
        boundary = 'bottom'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 2
    []
    [left_BC_4]
        type = FVPNMarshakRadiativeBC
        variable = phi_4
        boundary = 'left'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 4
    []
    [right_BC_4]
        type = FVPNMarshakRadiativeBC
        variable = phi_4
        boundary = 'right'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 4
    []
    [top_BC_4]
        type = FVPNMarshakRadiativeBC
        variable = phi_4
        boundary = 'top'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 4
    []
    [bottom_BC_4]
        type = FVPNMarshakRadiativeBC
        variable = phi_4
        boundary = 'bottom'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 4
    []
[]

[Executioner]
    type = Steady
[]

[Outputs]
    exodus = true
[]

