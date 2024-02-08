//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSFVOutflowPressureBC.h"

registerMooseObject("NavierStokesApp", INSFVOutflowPressureBC);

template <class T>
InputParameters
INSFVOutflowPressureBCTempl<T>::validParams()
{
  InputParameters params = FVDirichletBCBase::validParams();
  params += T::validParams();

  params.addClassDescription("Adds inlet boudnary conditon for the turbulent kinetic energy based "
                             "on turbulent intensity.");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "fluid density");
  params.addRequiredParam<MooseFunctorName>("mu", "Dynamic viscosity");
  return params;
}

template <class T>
INSFVOutflowPressureBCTempl<T>::INSFVOutflowPressureBCTempl(const InputParameters & params)
  : FVDirichletBCBase(params),
    T(params),
    _dim(_subproblem.mesh().dimension()),
    _u(getFunctor<ADReal>("u")),
    _v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _rho(getFunctor<ADReal>(NS::density)),
    _mu(getFunctor<ADReal>("mu"))
{
  if (_dim >= 2 && !_v)
    mooseError(
        "In two or more dimensions, the v velocity must be supplied using the 'v' parameter");
  if (_dim >= 3 && !_w)
    mooseError("In threedimensions, the w velocity must be supplied using the 'w' parameter");
}

template <class T>
ADReal
INSFVOutflowPressureBCTempl<T>::boundaryValue(const FaceInfo & fi) const
{
  // const auto boundary_face = singleSidedFaceArg(&fi);
  const auto state = determineState();

  const bool use_elem = fi.faceType(std::make_pair(_var.number(), _var.sys().number())) ==
                        FaceInfo::VarFaceNeighbors::ELEM;
  const auto elem_arg = use_elem ? makeElemArg(&fi.elem()) : makeElemArg(fi.neighborPtr());

  ADRealTensorValue grad_u;
  if (_dim == 1)
    grad_u = ADRealTensorValue(_u.gradient(elem_arg, state));
  else if (_dim == 2)
    grad_u =
        ADRealTensorValue(_u.gradient(elem_arg, state), (*_v).gradient(elem_arg, state));
  else
    grad_u = ADRealTensorValue(_u.gradient(elem_arg, state),
                               (*_v).gradient(elem_arg, state),
                               (*_w).gradient(elem_arg, state));

  const auto boundary_shear = (grad_u * fi.normal()) * fi.normal();

  return boundary_shear;
}
