//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVSP3ThermalRadiationDiffusion.h"
#include "MathUtils.h"

registerMooseObject("HeatTransferApp", FVSP3ThermalRadiationDiffusion);

InputParameters
FVSP3ThermalRadiationDiffusion::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params += FVDiffusionInterpolationInterface::validParams();
  params.addClassDescription(
      "Computes residual for the SP3 diffusion operator for finite volume method.");

  params.addRequiredParam<MooseFunctorName>("epsilon", "The optical thickness of the medium.");
  params.addRequiredParam<MooseFunctorName>("kappa", "The absorptivity of the medium.");

  MooseEnum order("first second", "first");
  params.addParam<MooseEnum>("order", order, "The order of the diffusion term.");

  MooseEnum coeff_interp_method("average harmonic", "harmonic");
  params.addParam<MooseEnum>(
      "coeff_interp_method",
      coeff_interp_method,
      "Switch that can select face interpolation method for diffusion coefficients.");

  params.set<unsigned short>("ghost_layers") = 2;

  return params;
}

FVSP3ThermalRadiationDiffusion::FVSP3ThermalRadiationDiffusion(const InputParameters & params)
  : FVFluxKernel(params),
    FVDiffusionInterpolationInterface(params),
    _optical_thickness(getFunctor<ADReal>("epsilon")),
    _absorptivity(getFunctor<ADReal>("kappa")),
    _order(getParam<MooseEnum>("order")),
    _coeff_interp_method(
      Moose::FV::selectInterpolationMethod(getParam<MooseEnum>("coeff_interp_method")))
{
    if (_order == "first")
      _mu_order = _mu_1_squared;
    else
      _mu_order = _mu_2_squared;
}

ADReal
FVSP3ThermalRadiationDiffusion::computeQpResidual()
{
  using namespace Moose::FV;
  const auto state = determineState();

  const auto dudn = gradUDotNormal(state, _correct_skewness);

  ADReal coef_face;

  // If we are on internal faces, we interpolate the diffusivity as usual
  if (_var.isInternalFace(*_face_info))
  {
    // Compute the diffusion coefficients at the element and the neighbor
    const auto optical_thickness_elem_squared = Utility::pow<2>(_optical_thickness(elemArg(), state));
    const auto optical_thickness_neighbor_squared = Utility::pow<2>(_optical_thickness(neighborArg(), state));
    const auto absorptivity_elem = (_absorptivity(elemArg(), state) > 1e-12 ? _absorptivity(elemArg(), state) : 1e-12);
    const auto absorptivity_neighbor = (_absorptivity(neighborArg(), state) > 1e-12 ? _absorptivity(neighborArg(), state) : 1e-12);;

    // Face interpolate coefficients
    interpolate(_coeff_interp_method,
                coef_face,
                optical_thickness_elem_squared / absorptivity_elem,
                optical_thickness_neighbor_squared / absorptivity_neighbor,
                *_face_info,
                true);
  }
  // Else we just use the boundary values (which depend on how the diffusion
  // coefficient is constructed)
  else
  {
    const auto face = singleSidedFaceArg();

    const auto optical_thickness_face_squared = Utility::pow<2>(_optical_thickness(face, state));
    const auto absorptivity_face = (_absorptivity(face, state) > 1e-12 ? _absorptivity(face, state) : 1e-12);

    coef_face = optical_thickness_face_squared / absorptivity_face;
  }

  // Scaling with order of the diffusion coefficient
  coef_face *= _mu_order;

  // Print for Debug
  // const auto facenorm = _face_info->normal();
  // const auto x_coord = _face_info->faceCentroid()(0);
  // const auto tempflux= -1 * coef_face * dudn;
  // if(coef_face != _mu_order) printf("ThermalDf %.3f: %.5f * %.5f = %.5f\n", x_coord, coef_face.value(), dudn.value(), tempflux.value());
  // int numorder = -1;
  // if(_order == "first") numorder = 1;
  // else if(_order == "second") numorder = 2;
  // if(x_coord >0.4 && x_coord < 0.47) printf("ThermalDf %.3f(%d) : %.5f\n", x_coord, numorder, tempflux.value());

  return -1 * coef_face * dudn;
}