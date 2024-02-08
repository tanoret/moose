//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PPhaseAdvection.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PPhaseAdvection);

InputParameters
WCNSFV2PPhaseAdvection::validParams()
{
  auto params = INSFVMassAdvection::validParams();
  return params;
}

WCNSFV2PPhaseAdvection::WCNSFV2PPhaseAdvection(const InputParameters & params)
  : INSFVMassAdvection(params)
{
}

ADReal
WCNSFV2PPhaseAdvection::computeQpResidual()
{
  const auto v = velocity();
  const auto rho_face = _rho(makeFace(*_face_info,
                                      limiterType(_advected_interp_method),
                                      MetaPhysicL::raw_value(v) * _normal > 0),
                             determineState());
  const auto var_face = _var(makeFace(*_face_info,
                                      limiterType(_advected_interp_method),
                                      MetaPhysicL::raw_value(v) * _normal > 0),
                             determineState());
  return _normal * v * rho_face * var_face;
}
