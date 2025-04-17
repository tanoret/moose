sigma_0 = 0.9
sigma_2 = 1.0
sigma_4 = 1.0
sigma_tr_1 = 1.0
sigma_tr_3 = 1.0
sigma_tr_5 = 1.0
source = 1.0

[Mesh]
    [cmg]
        type = CartesianMeshGenerator
        dim = 2
        dx = '5'
        dy = '5'
        ix = '100'
        iy = '100'
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
        expression = "if((x>=1.75 & x<=2.25 & y>=1.75 & y<=2.25)
                        |(x>=2.75 & x<=3.25 & y>=1.5 & y <=2.5)
                        |(x>=1.75 & x<=2.25 & y>=2.75 & y <=3.25)
                        |(x>=3.5 & x<=4.25 & y>=3.5 & y <=3.75)
                    , ${source}, 0)"
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

# [FVBCs]
#     [BC_0]
#         type = FVDirichletBC
#         boundary = 'left right bottom top'
#         variable = phi_0
#         value = 0
#     []
#     [BC_2]
#         type = FVDirichletBC
#         boundary = 'left right bottom top'
#         variable = phi_2
#         value = 0
#     []
#     [BC_4]
#         type = FVDirichletBC
#         boundary = 'left right bottom top'
#         variable = phi_4
#         value = 0
#     []
# []

[FVBCs]
    [BC_0]
        type = FVPNMarshakRadiativeBC
        variable = phi_0
        boundary = 'left right bottom top'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        n = 0
    []
    [BC_2]
        type = FVPNMarshakRadiativeBC
        variable = phi_2
        boundary = 'left right bottom top'
        phi_0 = 'phi_0'
        phi_2 = 'phi_2'
        phi_4 = 'phi_4'
        n = 2
    []
    [BC_4]
        type = FVPNMarshakRadiativeBC
        variable = phi_4
        boundary = 'left right bottom top'
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