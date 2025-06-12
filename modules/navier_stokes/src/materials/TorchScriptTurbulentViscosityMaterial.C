//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef LIBTORCH_ENABLED

#include "TorchScriptTurbulentViscosityMaterial.h"

registerMooseObject("MooseApp", TorchScriptTurbulentViscosityMaterial);

InputParameters
TorchScriptTurbulentViscosityMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Material object which relies on the evaluation of a TorchScript module to compute the turbulent dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<UserObjectName>(
      "torch_script_userobject",
      "The name of the user object which contains the torch script module.");
  params.addRequiredParam<MooseFunctorName>("k", "k value for Reynold's stress");
  params.addRequiredParam<MooseFunctorName>("eps", "epsilon for Reynold's stress");
  params.addParam<bool>("debug", "Value to control debugging console messages");
  params.addParam<Real>("mu_t_min", "Minimum value outputted for mu_t");

  return params;
}

TorchScriptTurbulentViscosityMaterial::TorchScriptTurbulentViscosityMaterial(const InputParameters & parameters)
  : Material(parameters),
    _mesh_dimension(_mesh.dimension()),
    _u_var(getFunctor<ADReal>("u")),
    _v_var(parameters.isParamValid("v") ? &(getFunctor<ADReal>("v")) : nullptr),
    _w_var(parameters.isParamValid("w") ? &(getFunctor<ADReal>("w")) : nullptr),
    _k(getFunctor<ADReal>("k")),
    _eps(getFunctor<ADReal>("eps")),
    _debug(getParam<bool>("debug")),
    _mu_t_min(getParam<Real>("mu_t_min")),
    _torch_script_userobject(getUserObject<TorchScriptUserObject>("torch_script_userobject")),
    _input_tensor(torch::zeros(
        {1, 2},
        torch::TensorOptions().dtype(torch::kFloat64).device(_app.getLibtorchDevice())))
    {
    _properties = &declareGenericPropertyByName<Real, false>("mu_t_torch");
}

void
TorchScriptTurbulentViscosityMaterial::initQpStatefulProperties()
{
  computeQpValues();
}

void
TorchScriptTurbulentViscosityMaterial::computeQpProperties()
{
  computeQpValues();
}

void
TorchScriptTurbulentViscosityMaterial::computeQpValues()
{
    
  const auto r = makeElemArg(_current_elem);
  const auto t = determineState();

  TensorValue<Real> grad_velocity;
  grad_velocity(0, 0) = MetaPhysicL::raw_value(_u_var.gradient(r,t)(0));
  if (_mesh_dimension > 1)
  {
    grad_velocity(0, 1) = MetaPhysicL::raw_value(_u_var.gradient(r,t)(1));
    grad_velocity(1, 0) = MetaPhysicL::raw_value(_v_var->gradient(r,t)(0));
    grad_velocity(1, 1) = MetaPhysicL::raw_value(_v_var->gradient(r,t)(1));
  }
  if (_mesh_dimension > 2)
  {
    grad_velocity(0, 2)= MetaPhysicL::raw_value(_u_var.gradient(r,t)(2));
    grad_velocity(1, 2)= MetaPhysicL::raw_value(_v_var->gradient(r,t)(2));
    grad_velocity(2, 0)= MetaPhysicL::raw_value(_w_var->gradient(r,t)(0));
    grad_velocity(2, 1)= MetaPhysicL::raw_value(_w_var->gradient(r,t)(1));
    grad_velocity(2, 2)= MetaPhysicL::raw_value(_w_var->gradient(r,t)(2));
  }

  const Real k = MetaPhysicL::raw_value(_k(r,t));
  const Real eps = std::max(MetaPhysicL::raw_value(_eps(r,t)), 1e-10);
  const Real timescale = (k / eps);

  const TensorValue<Real> sij = timescale * (grad_velocity + grad_velocity.transpose()) / 2.0;
  const TensorValue<Real> rij = timescale * (grad_velocity - grad_velocity.transpose()) / 2.0;

  const TensorValue<Real> I(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);

  const Real eta1 = std::min(std::max(sij.contract(sij), 1e-3), 1e3);
  const Real eta2 = std::min(std::max(rij.contract(rij), 1e-3), 1e3);

  const Real G1_term = eta1 ;
  const Real G2_term = sij.contract(sij * rij - rij * sij);
  const Real G3_term = sij.contract(sij * sij - 1./3. * I * I.contract(sij*sij));
  const Real k_term = 2./3. * k * sij.contract(I);

  auto input_accessor = _input_tensor.accessor<Real, 2>();
  input_accessor[0][0] =  eta1;
  input_accessor[0][1] =  eta2;

  const auto output = _torch_script_userobject.evaluate(_input_tensor);
  const auto output_accessor = output.accessor<Real, 2>();

  const auto G1 = output_accessor[0][0];
  const auto G2 = output_accessor[0][1];
  const auto G3 = output_accessor[0][2];

  const Real nu_t = (timescale / eta1) * (2. * k *(G1 * G1_term + G2 * G2_term + G3 * G3_term) + k_term);

  if (_debug)
  {
    _console << "-----------------------------------------" << std::endl;
    _console << "k: " << k << " eps: " << eps << std::endl;
    _console << "eta_1: " << input_accessor[0][0] << " eta_2: " << input_accessor[0][0] << std::endl;
    _console << "G1 term: " << G1_term << " G2 term: " << G2_term << " G3 term: " << G3_term << std::endl;
    _console << "G1: " << G1 << " G2: " << G2 << " G3: " << G3 << std::endl;
    _console << "nu_t output: " << nu_t << std::endl;
  }
  
  (*_properties)[_qp] = std::min(std::max(nu_t, 0.05), 0.15);
}

#endif
