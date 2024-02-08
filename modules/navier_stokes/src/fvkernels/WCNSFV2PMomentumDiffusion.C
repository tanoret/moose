//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMomentumDiffusion.h"
#include "INSFVRhieChowInterpolator.h"
#include "NS.h"
#include "NavierStokesMethods.h"
#include "SystemBase.h"
#include "NonlinearSystemBase.h"
#include "RelationshipManager.h"
#include "Factory.h"
#include "libmesh/nonlinear_solver.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMomentumDiffusion);

InputParameters
WCNSFV2PMomentumDiffusion::validParams()
{
  auto params = INSFVMomentumDiffusion::validParams();
  params.addRequiredParam<MooseFunctorName>("fd", "The phase fraction functor.");
  params.addClassDescription(
      "Implements the Laplace form of the two-phase viscous stress in the Navier-Stokes equation.");
  return params;
}

WCNSFV2PMomentumDiffusion::WCNSFV2PMomentumDiffusion(const InputParameters & params)
  : INSFVMomentumDiffusion(params),
    _alpha(getFunctor<ADReal>("fd"))
{
  if (!_u_var)
    paramError("u", "The u velocity must be defined for the divergence term in 'WCNSFV2PMomentumDiffusion'.");

  if (_dim >= 2 && !_v_var)
    paramError("v",
               "The v velocity must be defined for the divergence term in 'WCNSFV2PMomentumDiffusion'");

  if (_dim >= 3 && !_w_var)
    paramError("w",
               "The w velocity must be defined for the divergence term in 'WCNSFV2PMomentumDiffusion'");
}

ADReal
WCNSFV2PMomentumDiffusion::computeStrongResidual(const bool populate_a_coeffs)
{
  const Moose::StateArg state = determineState();
  const auto dudn = gradUDotNormal(state);
  ADReal face_mu;

  if (onBoundary(*_face_info))
    face_mu = _mu(makeCDFace(*_face_info), state);
  else
    Moose::FV::interpolate(_mu_interp_method,
                           face_mu,
                           _mu(elemArg(), state),
                           _mu(neighborArg(), state),
                           *_face_info,
                           true);

  // Protecting from negative viscosity at interpolation
  // to preserve convergence
  if (face_mu < 0.0)
  {
    if (!(_limit_interpolation))
      mooseWarning("Negative face viscosity has been encountered. Value ",
                   raw_value(face_mu),
                   " at ",
                   _face_info->faceCentroid(),
                   " limiting it to 0!");
    face_mu = 0;
  }

  if (populate_a_coeffs)
  {
    if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
        _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    {
      const auto dof_number = _face_info->elem().dof_number(_sys.number(), _var.number(), 0);
      _ae = dudn.derivatives()[dof_number];
      _ae *= -face_mu;
    }
    if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
        _face_type == FaceInfo::VarFaceNeighbors::BOTH)
    {
      const auto dof_number = _face_info->neighbor().dof_number(_sys.number(), _var.number(), 0);
      _an = dudn.derivatives()[dof_number];
      _an *= face_mu;
    }
  }

  Moose::FaceArg face;
  const bool skewness_correction =
      (_var.faceInterpolationMethod() == Moose::FV::InterpMethod::SkewCorrectedAverage);
  if (onBoundary(*_face_info))
    face = singleSidedFaceArg();
  else
    face = makeCDFace(*_face_info, skewness_correction);

  ADReal divergence = 0.0;
  if (_dim >= 1)
    divergence += _u_var->gradient(face, state)(0);
  if (_dim >= 2)
    divergence += _v_var->gradient(face, state)(1);
  if (_dim >= 3)
    divergence += _w_var->gradient(face, state)(2);

  ADReal dudn_transpose = 0.0;
  if (_complete_expansion)
  {
    ADRealTensorValue gradient;
    if (_dim == 1)
    {
      const auto & grad_u = _u_var->gradient(face, state);
      gradient = ADRealTensorValue(grad_u, ADRealVectorValue(0, 0, 0), ADRealVectorValue(0, 0, 0));
    }
    else if (_dim == 2)
    {
      const auto & grad_u = _u_var->gradient(face, state);
      const auto & grad_v = _v_var->gradient(face, state);
      gradient = ADRealTensorValue(grad_u, grad_v, ADRealVectorValue(0, 0, 0));
    }
    else // if (_dim == 3)
    {
      const auto & grad_u = _u_var->gradient(face, state);
      const auto & grad_v = _v_var->gradient(face, state);
      const auto & grad_w = _w_var->gradient(face, state);
      gradient = ADRealTensorValue(grad_u, grad_v, grad_w);
    }

    // Getting transpose of the gradient matrix
    const auto gradient_transpose = gradient.transpose();

    dudn_transpose += gradient_transpose.row(_index) * _face_info->normal();
  }

  return -face_mu * (dudn + dudn_transpose) + 2./3. * face_mu * divergence;
}
