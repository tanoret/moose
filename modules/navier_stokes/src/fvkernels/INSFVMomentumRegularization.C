//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumRegularization.h"
#include "INSFVRhieChowInterpolator.h"
#include "NS.h"
#include "SystemBase.h"

registerMooseObject("NavierStokesApp", INSFVMomentumRegularization);

InputParameters
INSFVMomentumRegularization::validParams()
{
  auto params = FVFluxKernel::validParams();
  params.addRequiredParam<MooseFunctorName>("filter_scaling", "alpha-regularization parameter for the NS equation");
  params.addClassDescription(
      "Implements a mesh-based regularization kernel for the Navier Stokes equation.");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVMomentumRegularization::INSFVMomentumRegularization(const InputParameters & params)
  : FVFluxKernel(params), _filter_scaling(getFunctor<ADReal>("filter_scaling"))
{
}

ADReal
INSFVMomentumRegularization::computeQpResidual()
{
  auto dudn = gradUDotNormal();

  Point _elem_centroid = _face_info->elemCentroid();
  Point _neighbor_centroid = _face_info->neighborCentroid();
  Real _cell_distance_norm = (_elem_centroid - _neighbor_centroid).norm();

  // Eventually, it will be nice to offer automatic-switching triggered by
  // input parameters to change between different interpolation methods for
  // this.
  const auto k = _filter_scaling(Moose::FV::makeCDFace(*_face_info, faceArgSubdomains()));

  return -1.0 * dudn * k * _cell_distance_norm;
}
