//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVSP3TemperatureBC.h"
#include "Function.h"
#include "HeatTransferModels.h"

registerMooseObject("HeatTransferApp", FVSP3TemperatureBC);

InputParameters
FVSP3TemperatureBC::validParams()
{
  InputParameters params = FVDirichletBCBase::validParams();
  params += FVDiffusionInterpolationInterface::validParams();
  params.addClassDescription("Semi-Transparent Boundary Condition for SP3 Radiation Temperature");

  params.addRequiredParam<MooseFunctorName>("Tb", "The temperature of the boundary");
  params.addParam<MooseFunctorName>("n1", 1.0, "The refraction coefficient of medium 1.");
  params.addParam<MooseFunctorName>("n2", 1.0, "The refraction coefficient of medium 2.");
  params.addRequiredParam<MooseFunctorName>("h", "The convective heat transfer coefficient of medium.");
  params.addRequiredParam<MooseFunctorName>("k", "The conductive heat transfer coefficient of medium.");
  params.addRequiredParam<MooseFunctorName>("epsilon","The optical thickness of the medium.");
  params.addRequiredParam<MooseFunctorName>("alpha", "The hemispheric emissivity of the medium.");
  params.addRequiredParam<Real>("nu1", "The maximum opaque frequency of the medium.");
  params.addRequiredParam<Real>( "nu_min", "The minimum frequency");

  return params;
}

FVSP3TemperatureBC::FVSP3TemperatureBC(const InputParameters & parameters)
  : FVDirichletBCBase(parameters),
    NeighborCoupleableMooseVariableDependencyIntermediateInterface(
        this, /*nodal=*/false, /*neighbor_nodal=*/false, /*is_fv=*/true),
    FVDiffusionInterpolationInterface(parameters),
    _var(*mooseVariableFV()),
    _Tb(getFunctor<ADReal>("Tb")),
    _n1(getFunctor<ADReal>("n1")),
    _n2(getFunctor<ADReal>("n2")),
    _h(getFunctor<ADReal>("h")),
    _k(getFunctor<ADReal>("k")),
    _epsilon(getFunctor<ADReal>("epsilon")),
    _alpha(getFunctor<ADReal>("alpha")),
    _nu1(getParam<Real>("nu1")),
    _nu_min(getParam<Real>("nu_min"))
{
  NeighborCoupleableMooseVariableDependencyIntermediateInterface::addMooseVariableDependency(&_var);

  if (_var.kind() == Moose::VarKindType::VAR_AUXILIARY)
    paramError("variable",
               "There should not be a need to specify a "
               "boundary condition for an auxiliary variable.");
}

ADReal
FVSP3TemperatureBC::boundaryValue(const FaceInfo & fi, const Moose::StateArg & state) const
{
  // Allow the functors to pick their side evaluation
  // const Moose::FaceArg face{
  //     &fi, Moose::FV::LimiterType::CentralDifference, true, false, nullptr, nullptr};
  auto face = singleSidedFaceArg(&fi);

  // Build the convective source at the boundary  
  const auto Tb = _Tb(face, state);
  const auto epsilon = _epsilon(face, state);
  const auto h = _h(face,  state);
  const auto k = _k(face,  state);
  const auto alpha = _alpha(face, state);
  printf("3\n");

  const auto T = _var(face, state);
  printf("4\n");
  const auto dudn = Moose::FV::gradUDotNormal(fi, _var, state, _correct_skewness);
  printf("5\n");

  const auto n1 = _n1(face, state);
  const Real abs_tol = 1E-8;
  const Real rel_tol = 1E-6;
  printf("6\n");

  const auto boundary_source = HeatTransferModels::integratedPlanckBand<ADReal>(n1, 1.0, Tb, _nu_min, _nu1, abs_tol, rel_tol); // put kappa = 1
  const auto cell_source = HeatTransferModels::integratedPlanckBand<ADReal>(n1, 1.0, T, _nu_min, _nu1, abs_tol, rel_tol); // put kappa = 1

  return  Tb - epsilon / h * k * dudn + alpha / h * Utility::pow<2>(_n2(face, state)/_n1(face, state)) * (boundary_source - cell_source) / (4.0);
}