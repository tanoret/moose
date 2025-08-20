bias_x = 0.7334168740991486
bias_y = 0.7334168740991486
bias_z = 0.7334168740991486

nx = 36
ny = 36
nz = 36
mid = 0.8
nmid = 30

[Mesh]
  [front_top_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_x = ${fparse bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse bias_z}
  []

  [front_top_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = ${fparse mid/2}
    zmax = 0.5
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
    xmax = ${fparse -mid/2}
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse bias_z}
  []

  [front_mid_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nmid}
    nz = ${nz}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_x = ${fparse bias_x}
    bias_z = ${fparse bias_z}
  []

  [front_mid_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${nmid}
    nz = ${nz}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_z = ${fparse bias_z}
  []

  [front_mid_left]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nmid}
    nz = ${nz}
    xmin = -0.5
    xmax = ${fparse -mid/2}
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_x = ${fparse 1/bias_x}
    bias_z = ${fparse bias_z}
  []

  [front_bottom_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_x = ${fparse bias_x}
    bias_y = ${fparse 1/bias_y}
    bias_z = ${fparse bias_z}
  []

  [front_bottom_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = ${fparse mid/2}
    zmax = 0.5
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
    xmax = ${fparse -mid/2}
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = ${fparse mid/2}
    zmax = 0.5
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse 1/bias_y}
    bias_z = ${fparse bias_z}
  []

  
  [mid_top_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nmid}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_x = ${fparse bias_x}
    bias_y = ${fparse bias_y}
  []

  [mid_top_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${ny}
    nz = ${nmid}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_y = ${fparse bias_y}
  []

  [mid_top_left]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nmid}
    xmin = -0.5
    xmax = ${fparse -mid/2}
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse bias_y}
  []

  [mid_mid_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nmid}
    nz = ${nmid}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_x = ${fparse bias_x}
  []

  [mid_mid_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${nmid}
    nz = ${nmid}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
  []

  [mid_mid_left]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nmid}
    nz = ${nmid}
    xmin = -0.5
    xmax = ${fparse -mid/2}
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_x = ${fparse 1/bias_x}

  []

  [mid_bottom_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nmid}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_x = ${fparse bias_x}
    bias_y = ${fparse 1/bias_y}
  []

  [mid_bottom_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${ny}
    nz = ${nmid}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_y = ${fparse 1/bias_y}
  []

  [mid_bottom_left]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nmid}
    xmin = -0.5
    xmax = ${fparse -mid/2}
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = ${fparse -mid/2}
    zmax = ${fparse mid/2}
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse 1/bias_y}
  []


  [back_top_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_x = ${fparse bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse 1/bias_z}
  []

  [back_top_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = -0.5
    zmax = ${fparse -mid/2}
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
    xmax = ${fparse -mid/2}
    ymin = ${fparse mid/2}
    ymax = 0.5
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse bias_y}
    bias_z = ${fparse 1/bias_z}
  []

  [back_mid_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nmid}
    nz = ${nz}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_x = ${fparse bias_x}
    bias_z = ${fparse 1/bias_z}
  []

  [back_mid_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${nmid}
    nz = ${nz}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_z = ${fparse 1/bias_z}
  []

  [back_mid_left]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${nmid}
    nz = ${nz}
    xmin = -0.5
    xmax = ${fparse -mid/2}
    ymin = ${fparse -mid/2}
    ymax = ${fparse mid/2}
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_x = ${fparse 1/bias_x}
    bias_z = ${fparse 1/bias_z}
  []

  [back_bottom_right]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nx}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse mid/2}
    xmax = 0.5
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_x = ${fparse bias_x}
    bias_y = ${fparse 1/bias_y}
    bias_z = ${fparse 1/bias_z}
  []

  [back_bottom_mid]
    type = GeneratedMeshGenerator
    dim = 3
    nx = ${nmid}
    ny = ${ny}
    nz = ${nz}
    xmin = ${fparse -mid/2}
    xmax = ${fparse mid/2}
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = -0.5
    zmax = ${fparse -mid/2}
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
    xmax = ${fparse -mid/2}
    ymin = -0.5
    ymax = ${fparse -mid/2}
    zmin = -0.5
    zmax = ${fparse -mid/2}
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse 1/bias_y}
    bias_z = ${fparse 1/bias_z}
  []
  
  
  [join_front_top]
      type = StitchedMeshGenerator
      inputs = 'front_top_right front_top_mid front_top_left'
      stitch_boundaries_pairs = 'left right; left right'
  []
  [join_front_mid]
    type = StitchedMeshGenerator
    inputs = 'front_mid_right front_mid_mid front_mid_left'
    stitch_boundaries_pairs = 'left right; left right'
  []
  [join_front_bottom]
    type = StitchedMeshGenerator
    inputs = 'front_bottom_right front_bottom_mid front_bottom_left'
    stitch_boundaries_pairs = 'left right; left right'
  []

  [join_mid_top]
    type = StitchedMeshGenerator
    inputs = 'mid_top_right mid_top_mid mid_top_left'
    stitch_boundaries_pairs = 'left right; left right'
  []
  [join_mid_mid]
    type = StitchedMeshGenerator
    inputs = 'mid_mid_right mid_mid_mid mid_mid_left'
    stitch_boundaries_pairs = 'left right; left right'
  []
  [join_mid_bottom]
    type = StitchedMeshGenerator
    inputs = 'mid_bottom_right mid_bottom_mid mid_bottom_left'
    stitch_boundaries_pairs = 'left right; left right'
  []

  [join_back_top]
    type = StitchedMeshGenerator
    inputs = 'back_top_right back_top_mid back_top_left'
    stitch_boundaries_pairs = 'left right; left right'
  []
  [join_back_mid]
    type = StitchedMeshGenerator
    inputs = 'back_mid_right back_mid_mid back_mid_left'
    stitch_boundaries_pairs = 'left right; left right'
  []
  [join_back_bottom]
    type = StitchedMeshGenerator
    inputs = 'back_bottom_right back_bottom_mid back_bottom_left'
    stitch_boundaries_pairs = 'left right; left right'
  []

  
  [join_front]
    type = StitchedMeshGenerator
    inputs = 'join_front_top join_front_mid join_front_bottom'
    stitch_boundaries_pairs = 'bottom top; bottom top'
  []
  [join_mid]
    type = StitchedMeshGenerator
    inputs = 'join_mid_top join_mid_mid join_mid_bottom'
    stitch_boundaries_pairs = 'bottom top; bottom top'
  []
  [join_back]
    type = StitchedMeshGenerator
    inputs = 'join_back_top join_back_mid join_back_bottom'
    stitch_boundaries_pairs = 'bottom top; bottom top'
  []

  [joint_all]
    type = StitchedMeshGenerator
    inputs = 'join_front join_mid join_back'
    stitch_boundaries_pairs = 'back front; back front'
  []
[]

# [Mesh]
#   [front_top_rignt]
#     type = GeneratedMeshGenerator
#     dim = 3
#     nx = ${nx}
#     ny = ${ny}
#     nz = ${nz}
#     xmin = 0.0
#     xmax = 0.5
#     ymin = 0.0
#     ymax = 0.5
#     zmin = 0.0
#     zmax = 0.5
#     bias_x = ${fparse bias_x}
#     bias_y = ${fparse bias_y}
#     bias_z = ${fparse bias_z}
#   []

#   [front_top_left]
#     type = GeneratedMeshGenerator
#     dim = 3
#     nx = ${nx}
#     ny = ${ny}
#     nz = ${nz}
#     xmin = -0.5
#     xmax = 0.0
#     ymin = 0.0
#     ymax = 0.5
#     zmin = 0.0
#     zmax = 0.5
#     bias_x = ${fparse 1/bias_x}
#     bias_y = ${fparse bias_y}
#     bias_z = ${fparse bias_z}
#   []

#   [front_bottom_right]
#     type = GeneratedMeshGenerator
#     dim = 3
#     nx = ${nx}
#     ny = ${ny}
#     nz = ${nz}
#     xmin = 0.0
#     xmax = 0.5
#     ymin = -0.5
#     ymax = 0.0
#     zmin = 0.0
#     zmax = 0.5
#     bias_x = ${fparse bias_x}
#     bias_y = ${fparse 1/bias_y}
#     bias_z = ${fparse bias_z}
#   []
#   [front_bottom_left]
#       type = GeneratedMeshGenerator
#       dim = 3
#       nx = ${nx}
#       ny = ${ny}
#       nz = ${nz}
#       xmin = -0.5
#       xmax = 0.0
#       ymin = -0.5
#       ymax = 0.0
#       zmin = 0.0
#       zmax = 0.5
#       bias_x = ${fparse 1/bias_x}
#       bias_y = ${fparse 1/bias_y}
#       bias_z = ${fparse bias_z}
#   []
#   [back_top_right]
#     type = GeneratedMeshGenerator
#     dim = 3
#     nx = ${nx}
#     ny = ${ny}
#     nz = ${nz}
#     xmin = 0.0
#     xmax = 0.5
#     ymin = 0.0
#     ymax = 0.5
#     zmin = -0.5
#     zmax = 0.0
#     bias_x = ${fparse bias_x}
#     bias_y = ${fparse bias_y}
#     bias_z = ${fparse 1/bias_z}
#   []
#   [back_top_left]
#       type = GeneratedMeshGenerator
#       dim = 3
#       nx = ${nx}
#       ny = ${ny}
#       nz = ${nz}
#       xmin = -0.5
#       xmax = 0.0
#       ymin = 0.0
#       ymax = 0.5
#       zmin = -0.5
#       zmax = 0.0
#       bias_x = ${fparse 1/bias_x}
#       bias_y = ${fparse bias_y}
#       bias_z = ${fparse 1/bias_z}
#   []
#   [back_bottom_right]
#       type = GeneratedMeshGenerator
#       dim = 3
#       nx = ${nx}
#       ny = ${ny}
#       nz = ${nz}
#       xmin = 0.0
#       xmax = 0.5
#       ymin = -0.5
#       ymax = 0.0
#       zmin = -0.5
#       zmax = 0.0
#       bias_x = ${fparse bias_x}
#       bias_y = ${fparse 1/bias_y}
#       bias_z = ${fparse 1/bias_z}
#   []
#   [back_bottom_left]
#       type = GeneratedMeshGenerator
#       dim = 3
#       nx = ${nx}
#       ny = ${ny}
#       nz = ${nz}
#       xmin = -0.5
#       xmax = 0.0
#       ymin = -0.5
#       ymax = 0.0
#       zmin = -0.5
#       zmax = 0.0
#       bias_x = ${fparse 1/bias_x}
#       bias_y = ${fparse 1/bias_y}
#       bias_z = ${fparse 1/bias_z}
#   []
#   [join_front_top]
#       type = StitchedMeshGenerator
#       inputs = 'front_top_rignt front_top_left'
#       stitch_boundaries_pairs = 'left right'
#   []
#   [join_front_bottom]
#     type = StitchedMeshGenerator
#     inputs = 'front_bottom_right front_bottom_left'
#     stitch_boundaries_pairs = 'left right'
#   []
#   [join_back_top]
#     type = StitchedMeshGenerator
#     inputs = 'back_top_right back_top_left'
#     stitch_boundaries_pairs = 'left right'
#   []
#   [join_back_bottom]
#     type = StitchedMeshGenerator
#     inputs = 'back_bottom_right back_bottom_left'
#     stitch_boundaries_pairs = 'left right'
#   []
#   [join_front]
#     type = StitchedMeshGenerator
#     inputs = 'join_front_top join_front_bottom'
#     stitch_boundaries_pairs = 'bottom top'
#   []
#   [join_back]
#     type = StitchedMeshGenerator
#     inputs = 'join_back_top join_back_bottom'
#     stitch_boundaries_pairs = 'bottom top'
#   []
#   [joint_all]
#     type = StitchedMeshGenerator
#     inputs = 'join_front join_back'
#     stitch_boundaries_pairs = 'back front'
#   []
# []