//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVPNThermalRadiation.h"
#include "MathUtils.h"
#include "HeatConductionNames.h"

registerMooseObject("HeatTransferApp", FVPNThermalRadiation);

InputParameters
FVPNThermalRadiation::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params += FVDiffusionInterpolationInterface::validParams();
  params.addClassDescription(
      "Kernel that assembles the equations for PN thermal radiation transport.");

  params.addRequiredParam<unsigned int>("n", "N-th Order of the Equation.");

  params.addParam<MooseFunctorName>(
      "phi_n_minus_2", 0.0, "Lagged order of the radiation flux.");
  params.addParam<MooseFunctorName>(
      "phi_n_plus_2", 0.0, "Leading order of the radiation flux.");

  params.addParam<MooseFunctorName>(
        "sigma_n_minus_1", 1.0, "Lagged order of the radiation flux.");
  params.addParam<MooseFunctorName>(
        "sigma_n_plus_1", 1.0, "Leading order of the radiation flux.");

  MooseEnum coeff_interp_method("average harmonic", "harmonic");
  params.addParam<MooseEnum>(
      "coeff_interp_method",
      coeff_interp_method,
      "Switch that can select face interpolation method for diffusion coefficients.");

  return params;
}

FVPNThermalRadiation::FVPNThermalRadiation(const InputParameters & params)
  : FVFluxKernel(params),
    FVDiffusionInterpolationInterface(params),
    _n(getParam<unsigned int>("n")),
    _phi_n_minus_2(getFunctor<ADReal>("phi_n_minus_2")),
    _phi_n_plus_2(getFunctor<ADReal>("phi_n_plus_2")),
    _sigma_n_minus_1(getFunctor<ADReal>("sigma_n_minus_1")),
    _sigma_n_plus_1(getFunctor<ADReal>("sigma_n_plus_1"))
{
  const auto & interp_method = getParam<MooseEnum>("coeff_interp_method");
  if (interp_method == "average")
    _coeff_interp_method = Moose::FV::InterpMethod::Average;
  else if (interp_method == "harmonic")
    _coeff_interp_method = Moose::FV::InterpMethod::HarmonicAverage;
}

ADReal
FVPNThermalRadiation::computeQpResidual()
{
  using namespace Moose::FV;
  const auto state = determineState();
  const auto elem_arg = elemArg();
  const auto rn = Real(_n);

  // Building the residual for the leading order
  const auto n_order_plus_2 = (rn+1)*(rn+2)/((2*rn+1)*(2*rn+3));
  ADReal sigma_n_plus_1_interpolated;
  interpolate(_coeff_interp_method,
    sigma_n_plus_1_interpolated,
    _sigma_n_plus_1(elemArg(), state),
    _sigma_n_plus_1(neighborArg(), state),
    *_face_info,
    true);
  const auto grad_phi_n_plus_2_normal = _phi_n_plus_2.gradient(elem_arg, state) * _face_info->normal();
  const auto res_leading_order = -n_order_plus_2 / sigma_n_plus_1_interpolated * grad_phi_n_plus_2_normal;

  // Building the residual for the leagged order
  const auto n_order_minus_2 = rn*(rn-1)/((2*rn+1)*(2*rn-1));
  ADReal sigma_n_minus_1_interpolated;
  interpolate(_coeff_interp_method,
    sigma_n_minus_1_interpolated,
    _sigma_n_minus_1(elemArg(), state),
    _sigma_n_minus_1(neighborArg(), state),
    *_face_info,
    true);
  const auto grad_phi_n_minus_2_normal = _phi_n_minus_2.gradient(elem_arg, state) * _face_info->normal();
  const auto res_lagged_order = -n_order_minus_2 / sigma_n_minus_1_interpolated * grad_phi_n_minus_2_normal;

  // Building the residual for the current term
  const auto current_n_order_plus_2 = Utility::pow<2>(rn+1)/((2*rn+1)*(2*rn+3));
  const auto current_n_order_minus_2 = Utility::pow<2>(rn)/((2*rn+1)*(2*rn-1));
  const auto res_order = -(current_n_order_plus_2/sigma_n_plus_1_interpolated + 
                           current_n_order_minus_2/sigma_n_minus_1_interpolated) * gradUDotNormal(state, _correct_skewness);

  return res_leading_order + res_lagged_order + res_order;
}
