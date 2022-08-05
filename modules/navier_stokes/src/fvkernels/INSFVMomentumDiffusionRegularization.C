//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVMomentumDiffusionRegularization.h"
#include "INSFVRhieChowInterpolator.h"
#include "NS.h"
#include "SystemBase.h"
#include "RelationshipManager.h"
#include "Factory.h"

registerMooseObject("NavierStokesApp", INSFVMomentumDiffusionRegularization);

InputParameters
INSFVMomentumDiffusionRegularization::validParams()
{
  auto params = INSFVFluxKernel::validParams();
  params.addRequiredParam<MooseFunctorName>("alpha", "The regularization parameter");
  params.addClassDescription(
      "Implements the Laplace form of the viscous stress in the Navier-Stokes equation.");
  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

INSFVMomentumDiffusionRegularization::INSFVMomentumDiffusionRegularization(const InputParameters & params)
  : FVFluxKernel(params), _alpha(getFunctor<ADReal>("alpha"))
{
  // We have no specific releationship managers for regularization diffusion kernels.
}

ADReal
INSFVMomentumDiffusionRegularization::computeQpResidual()
{
  // normal gradient of velocity
  const auto dudn = gradUDotNormal();

  // regularization factor at the face which element is being assembled
  const auto face_alpha = _alpha(Moose::FV::makeCDFace(*_face_info, faceArgSubdomains()));

  // norm of the ditance between the cell that is being assembled and its neighbour at the face
  Real _cell_distance_norm = (_face_info->elemCentroid() - _face_info->neighborCentroid()).norm();

  return -face_alpha * _cell_distance_norm * dudn;
}
