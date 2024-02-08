//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PLevelSet.h"
#include "NS.h"
#include "MathFVUtils.h"
#include "NavierStokesMethods.h"

registerMooseObject("NavierStokesApp", WCNSFV2PLevelSet);

InputParameters
WCNSFV2PLevelSet::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription("Computes the level set source term.");
  params.addParam<MooseFunctorName>("interface_thickness", 1e-3, "The interface thickness.");
  params.addParam<MooseFunctorName>("gamma", 1.0, "The reinitialization parameter.");
  params.addParam<MooseFunctorName>("regularizer_control", 1.0, "The reinitialization parameter.");

  params.set<unsigned short>("ghost_layers") = 2;
  return params;
}

WCNSFV2PLevelSet::WCNSFV2PLevelSet(const InputParameters & params)
  : FVFluxKernel(params),
    _reinitialization_parameter(getFunctor<ADReal>("gamma")),
    _interface_thickness(getFunctor<ADReal>("interface_thickness")),
    _regularizer_control(getFunctor<ADReal>("regularizer_control"))
{
}

ADReal
WCNSFV2PLevelSet::computeQpResidual()
{
  constexpr Real offset = 1e-15; // prevents explosion of sqrt(x) derivative to infinity

  const auto face = onBoundary(*_face_info) ? singleSidedFaceArg() : makeCDFace(*_face_info);
  const auto state = determineState();
  // const auto old_state = Moose::StateArg(1, Moose::SolutionIterationType::Nonlinear);

  const auto phi = _var(face, state);
  const auto grad_phi = _var.gradient(face, state);
  const auto grad_phi_norm = NS::computeSpeed(_var.gradient(face, state)) + offset;

  const auto regularizer =
      _regularizer_control(face, state) * phi * (1.0 - phi) * gradUDotNormal(state) / grad_phi_norm;

  const auto dPhidn = _interface_thickness(face, state) * gradUDotNormal(state);

  // ADReal rho_cp_face;
  // if (onBoundary(*_face_info))
  // {
  //   const auto ssf = singleSidedFaceArg();
  //   rho_cp_face = _rho(ssf, state) * _cp(ssf, state);
  // }
  // else
  // {
  //   // Interpolate the heat capacity
  //   const auto face_elem = elemArg();
  //   const auto face_neighbor = neighborArg();
  //   interpolate(Moose::FV::InterpMethod::Average,
  //               rho_cp_face,
  //               _rho(face_elem, state) * _cp(face_elem, state),
  //               _rho(face_neighbor, state) * _cp(face_neighbor, state),
  //               *_face_info,
  //               true);
  // }

  return -1 * _reinitialization_parameter(face, state) * (dPhidn - regularizer);
}
