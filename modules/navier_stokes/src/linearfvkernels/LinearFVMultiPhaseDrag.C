//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMultiPhaseDrag.h"
#include "Assembly.h"
#include "SubProblem.h"

registerMooseObject("MooseApp", LinearFVMultiPhaseDrag);

InputParameters
LinearFVMultiPhaseDrag::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription(
      "Represents the matrix and right hand side contributions of a reaction "
      "term ($c u$) in a partial differential equation.");

  params.addRequiredParam<MooseFunctorName>("drag_coeff", "The drag coefficient.");
  params.addRequiredParam<MooseFunctorName>("coupled_velocity", "The coupled velocity for drag.");

  return params;
}

LinearFVMultiPhaseDrag::LinearFVMultiPhaseDrag(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _drag_coefficient(getFunctor<Real>("drag_coeff")),
    _coupled_velocity(getFunctor<Real>("coupled_velocity"))
{
}

Real
LinearFVMultiPhaseDrag::computeMatrixContribution()
{
  return _drag_coefficient(makeElemArg(_current_elem_info->elem()), determineState()) *
         _current_elem_volume;
}

Real
LinearFVMultiPhaseDrag::computeRightHandSideContribution()
{
  return _drag_coefficient(makeElemArg(_current_elem_info->elem()), determineState()) *
         _coupled_velocity(makeElemArg(_current_elem_info->elem()), determineState()) *
         _current_elem_volume;
}
