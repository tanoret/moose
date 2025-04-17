//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVSP3ThermalRadiationSourceSink.h"
#include "MathUtils.h"
#include "HeatConductionNames.h"

registerMooseObject("HeatTransferApp", FVSP3ThermalRadiationSourceSink);

InputParameters
FVSP3ThermalRadiationSourceSink::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params.addClassDescription("Elemental kernel to the thermal radiation source and sink.");

  params.addRequiredRangeCheckedParam<MooseFunctorName>("T", "T>0", "The temperature of the medium.");
  params.addRequiredRangeCheckedParam<MooseFunctorName>("nu", "T>0", "The mean frequency of the thermal radiation band.");
  params.addParam<MooseFunctorName>("refraction_index", 1.0, "The refraction index in the spectral band.");
  params.addRequiredParam<MooseFunctorName>("kappa", "The absorptivity of the medium.");

  return params;
}

FVSP3ThermalRadiationSourceSink::FVSP3ThermalRadiationSourceSink(const InputParameters & params)
  : FVElementalKernel(params),
    _T(getFunctor<ADReal>("T")),
    _nu(getFunctor<ADReal>("nu")),
    _n1(getFunctor<ADReal>("refraction_index")),
    _absorptivity(getFunctor<ADReal>("kappa"))
{
}

ADReal
FVSP3ThermalRadiationSourceSink::computeQpResidual()
{
  // Convenient arguments
  const auto state = determineState();
  const auto elem_arg = makeElemArg(_current_elem);

  // Build the emission source
  const auto n1_pow_2 = Utility::pow<2>(_n1(elem_arg, state));
  const auto nu = _nu(elem_arg, state);
  const auto nu_pow_3 = Utility::pow<3>(nu);
  const auto T = _T(elem_arg, state);

  const auto pre_factor = n1_pow_2 * 2.0 * HeatConduction::Constants::hp * nu_pow_3 / (Utility::pow<2>(HeatConduction::Constants::c));
  const auto inv_thermal_source = std::exp(HeatConduction::Constants::hp*nu/(HeatConduction::Constants::kb * T)) - 1.0;
  const auto thermal_rad_source = 4.0 * libMesh::pi * _absorptivity(elem_arg, state) * pre_factor / inv_thermal_source;

  // Build the absorption sink
  const auto thermal_rad_sink = _absorptivity(elem_arg, state) * _var(elem_arg, state);

  // Return the residual
  return thermal_rad_sink - thermal_rad_source;
}
