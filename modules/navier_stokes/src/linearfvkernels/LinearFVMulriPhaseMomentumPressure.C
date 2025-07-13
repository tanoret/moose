//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMultiPhaseMomentumPressure.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "NS.h"
#include "FEProblemBase.h"

registerMooseObject("NavierStokesApp", LinearFVMultiPhaseMomentumPressure);

InputParameters
LinearFVMultiPhaseMomentumPressure::validParams()
{
  InputParameters params = LinearFVMomentumPressure::validParams();
  params.addClassDescription("Represents the pressure gradient term in the Navier Stokes momentum "
                             "equations in two phases, added to the right hand side.");
  params.addRequiredParam<MooseFunctorName>("alpha", "The phase fraction.");
  return params;
}

LinearFVMultiPhaseMomentumPressure::LinearFVMultiPhaseMomentumPressure(const InputParameters & params)
  : LinearFVMomentumPressure(params),
    _alpha(getFunctor<Real>("alpha"))
{
}

Real
LinearFVMultiPhaseMomentumPressure::computeRightHandSideContribution()
{
  const auto alpha = _alpha(makeElemArg(_current_elem_info->elem()), determineState());
  const auto dof_value = _current_elem_info->dofIndices()[_pressure_sys_num][_pressure_var_num];
  return -alpha*(*_pressure_gradient[_index])(dof_value)*_current_elem_volume;
}
