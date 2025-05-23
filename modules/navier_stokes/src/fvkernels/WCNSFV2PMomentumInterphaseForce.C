//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PMomentumInterphaseForce.h"

#include "NS.h"

registerMooseObject("NavierStokesApp", WCNSFV2PMomentumInterphaseForce);

InputParameters
WCNSFV2PMomentumInterphaseForce::validParams()
{
  InputParameters params = FVElementalKernel::validParams();
  params += INSFVMomentumResidualObject::validParams();

  params.addClassDescription("Computes the interphase drag force.");

  params.addRequiredParam<MooseFunctorName>("u_var_coupled", "The velocity of the coupled phase in the x direction.");
  params.addRequiredParam<MooseFunctorName>("v_var_coupled", "The velocity of the coupled phase in the y direction.");
  params.addRequiredParam<MooseFunctorName>("w_var_coupled", "The velocity of the coupled phase in the z direction.");


  params.addRequiredParam<MooseFunctorName>("u_var", "The velocity of the phase in the x direction.");
  params.addRequiredParam<MooseFunctorName>("v_var", "The velocity of the phase in the y direction.");
  params.addRequiredParam<MooseFunctorName>("w_var", "The velocity of the phase in the z direction.");

  params.addRequiredParam<MooseFunctorName>("fd", "Dispersed phase fraction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "Phase density.");

  params.addParam<bool>("activate_drag", true, "Activate drag force.");
  params.addParam<bool>("activate_lift", false, "Activate lift force.");
  params.addParam<bool>("activate_virtual_mass", false, "Activate virtual mass force.");

  params.addRequiredParam<MooseFunctorName>("drag_coef", "The drag coefficient.");

  return params;
}

WCNSFV2PMomentumInterphaseForce::WCNSFV2PMomentumInterphaseForce(const InputParameters & params)
  : FVElementalKernel(params),
    INSFVMomentumResidualObject(*this),
    _dim(_subproblem.mesh().dimension()),
    _u_var_coupled(params.isParamValid("u_var_coupled")
                   ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("u_var_coupled", 0))
                   : nullptr),
    _v_var_coupled(params.isParamValid("v_var_coupled")
                   ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("v_var_coupled", 0))
                   : nullptr),
    _w_var_coupled(params.isParamValid("w_var_coupled")
                   ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("w_var_coupled", 0))
                   : nullptr),
    _u_var(params.isParamValid("u_var")
                   ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("u_var", 0))
                   : nullptr),
    _v_var(params.isParamValid("v_var")
                   ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("v_var", 0))
                   : nullptr),
    _w_var(params.isParamValid("w_var")
                   ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("w_var", 0))
                   : nullptr),
    _drag_coef(getFunctor<ADReal>("drag_coef")),
    _alpha(getFunctor<ADReal>("fd")),
    _rho(getFunctor<ADReal>(NS::density)),
    _bool_activate_drag(getParam<bool>("activate_drag")),
    _bool_activate_lift(getParam<bool>("activate_lift")),
    _bool_activate_virtual_mass(getParam<bool>("activate_virtual_mass"))
{
  if (_index == 0 && !_u_var_coupled)
    paramError("u_var_coupled", "The coupled veloicty in the 'x' direction must be provided for 'x' momentum component");

  if (_index == 1 && !_v_var_coupled)
    paramError("v_var_coupled", "The coupled veloicty in the 'y' direction must be provided for 'y' momentum component");

  if (_index == 2 && !_w_var_coupled)
    paramError("z_var_coupled", "The coupled veloicty in the 'z' direction must be provided for 'z' momentum component");

  if (_bool_activate_lift)
  {
    if (!_u_var)
      paramError("u_var", "The pahse velocity in the 'x' direction must be provided if the lift force is active.");
    
    if (_dim > 1 && !_v_var)
      paramError("v_var", "The pahse velocity in the 'y' direction must be provided if the lift force is active.");

    if (_dim > 2 && !_w_var)
      paramError("w_var", "The pahse velocity in the 'z' direction must be provided if the lift force is active.");
  }
}

ADReal
WCNSFV2PMomentumInterphaseForce::computeQpResidual()
{
  const auto r = makeElemArg(_current_elem);
  const auto t = determineState();

  ADReal force = 0.;

  if (_bool_activate_drag)
  {
    ADReal dragvelocity;
    switch (_index)
    {
      case 0:
        dragvelocity = (*_u_var_coupled)(r,t);
        break;
      case 1:
        dragvelocity = (*_v_var_coupled)(r,t);
        break;
      case 2:
        dragvelocity = (*_w_var_coupled)(r,t);
        break;
    }

    force += _drag_coef(r, t) * (dragvelocity - _var(r, t));
  }

  if (_bool_activate_lift)
  {

    ADRealVectorValue curl = 0;
    if (_dim > 1)
    {
      const auto u_grad = _u_var->gradient(r, t);
      const auto v_grad = _v_var->gradient(r, t);
      curl(2) = v_grad(0) - u_grad(1);
      if (_dim > 2)
      {
        const auto w_grad = _w_var->gradient(r, t);
        curl(0) = w_grad(1) - v_grad(2);
        curl(1) = u_grad(2) - w_grad(0);
      }

      ADRealVectorValue speed_diff = 0;
      speed_diff(0) = (*_u_var_coupled)(r, t) - (*_u_var)(r, t);
      if (_dim > 1)
      {
          speed_diff(1) = (*_v_var_coupled)(r, t) - (*_v_var)(r, t);
          if (_dim > 2)
              speed_diff(2) = (*_w_var_coupled)(r, t) - (*_w_var)(r, t);
      }

      const auto cross_term = (speed_diff.cross(curl))(_index);
      force += 0.5 * _alpha(r, t) * _rho(r, t) * cross_term;
    }
  }

  if (_bool_activate_virtual_mass)
  {
    // Adding transient term
    ADRealVectorValue term_advection_phase = 0;
    ADRealVectorValue term_transient_phase = 0;
    ADRealVectorValue term_advection_coupled_phase = 0;
    ADRealVectorValue term_transient_coupled_phase = 0;

    std::vector<ADRealVectorValue *> term_transient = {&term_transient_phase, &term_transient_coupled_phase};
    std::vector<ADRealVectorValue *> term_advection = {&term_advection_phase, &term_advection_coupled_phase};
    std::vector<const Moose::Functor<ADReal> *> vel_phase = {_u_var, _v_var, _w_var};
    std::vector<const Moose::Functor<ADReal> *> vel_phase_coupled = {_u_var_coupled, _v_var_coupled, _w_var_coupled};
    std::vector<std::vector<const Moose::Functor<ADReal> *>>vel = {vel_phase, vel_phase_coupled};

    for (std::size_t i = 0; i < 2; ++i)
    {
      if (_subproblem.isTransient())
      {
        (*term_transient[i])(0) += vel[i][0]->dot(r, t);
        if (_dim > 1)
          (*term_transient[i])(1) += vel[i][1]->dot(r, t);
        if (_dim > 2)
          (*term_transient[i])(2) += vel[i][2]->dot(r, t);
      }

      // Adding advection term
      const auto u_velocity = (*vel[i][0])(r, t);
      const auto u_grad = vel[i][0]->gradient(r, t);
      (*term_advection[i])(0) += u_velocity * u_grad(0);
      if (_dim > 1)
      {
        const auto v_velocity = (*vel[i][1])(r, t);
        const auto v_grad = vel[i][1]->gradient(r, t);
        (*term_advection[i])(0) += v_velocity * u_grad(1);
        (*term_advection[i])(1) += u_velocity * v_grad(0) + v_velocity * v_grad(1);
        if (_dim > 2)
        {
          const auto w_velocity = (*vel[i][2])(r, t);
          const auto w_grad = vel[i][2]->gradient(r, t);
          (*term_advection[i])(0) += w_velocity * u_grad(2);
          (*term_advection[i])(1) += w_velocity * v_grad(2);
          (*term_advection[i])(2) +=
              u_velocity * w_grad(0) + v_velocity * w_grad(1) + w_velocity * w_grad(2);
        }
      }
    }

    ADReal relative_acceleration = (*term_advection[0])(_index) - (*term_advection[1])(_index);
    if (_subproblem.isTransient())
      relative_acceleration += (*term_transient[0])(_index) - (*term_transient[1])(_index);

    force += 0.5 * _alpha(r, t) * _rho(r, t) * relative_acceleration;
  }

  return force;
}
