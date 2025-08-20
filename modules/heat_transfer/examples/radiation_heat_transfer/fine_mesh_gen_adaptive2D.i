bias_x = 0.7334168740991486
bias_y = 0.7334168740991486
bias_z = 0.7334168740991486

nx = 26
ny = 26
nz = 26
mid = 0.8
nmid = 30

[Mesh]
  [front_top_rignt]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmin = 0.0
    xmax = 0.5
    ymin = 0.0
    ymax = 0.5
    bias_x = ${fparse bias_x}
    bias_y = ${fparse bias_y}
  []

  [front_top_left]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmin = -0.5
    xmax = 0.0
    ymin = 0.0
    ymax = 0.5
    bias_x = ${fparse 1/bias_x}
    bias_y = ${fparse bias_y}
  []

  [front_bottom_right]
    type = GeneratedMeshGenerator
    dim = 2
    nx = ${nx}
    ny = ${ny}
    xmin = 0.0
    xmax = 0.5
    ymin = -0.5
    ymax = 0.0
    bias_x = ${fparse bias_x}
    bias_y = ${fparse 1/bias_y}
  []
  [front_bottom_left]
      type = GeneratedMeshGenerator
      dim = 2
      nx = ${nx}
      ny = ${ny}
      xmin = -0.5
      xmax = 0.0
      ymin = -0.5
      ymax = 0.0
      bias_x = ${fparse 1/bias_x}
      bias_y = ${fparse 1/bias_y}
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
 
  [join_front]
    type = StitchedMeshGenerator
    inputs = 'join_front_top join_front_bottom'
    stitch_boundaries_pairs = 'bottom top'
  []
[]