//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifdef LIBTORCH_ENABLED

#include "TorchScriptTurbulentAnisotropyMaterial.h"
#include "NS.h"

registerMooseObject("MooseApp", TorchScriptTurbulentAnisotropyMaterial);

InputParameters
TorchScriptTurbulentAnisotropyMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Material object which relies on the evaluation of a TorchScript module to compute the turbulent dynamic viscosity.");
  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "The density.");
  params.addRequiredParam<UserObjectName>(
      "torch_script_userobject",
      "The name of the user object which contains the torch script module.");
  params.addRequiredParam<MooseFunctorName>("k", "k value for Reynold's stress");
  params.addRequiredParam<MooseFunctorName>("eps", "epsilon for Reynold's stress");
  params.addRequiredParam<bool>("debug", "Value to control debugging console messages");
  params.addRequiredParam<bool>("use_NN", "Value to control whether using NN or ASRM");
  params.addParam<MooseFunctorName>("property_prefix", "b", "Anisotropy name prefix ex. b");

  return params;
}

TorchScriptTurbulentAnisotropyMaterial::TorchScriptTurbulentAnisotropyMaterial(const InputParameters & parameters)
  : Material(parameters),
    _mesh_dimension(_mesh.dimension()),
    _u_var(getFunctor<ADReal>("u")),
    _v_var(parameters.isParamValid("v") ? &(getFunctor<ADReal>("v")) : nullptr),
    _w_var(parameters.isParamValid("w") ? &(getFunctor<ADReal>("w")) : nullptr),
    _rho(getFunctor<ADReal>(NS::density)),
    _k(getFunctor<ADReal>("k")),
    _eps(getFunctor<ADReal>("eps")),
    _debug(getParam<bool>("debug")),
    _use_NN(getParam<bool>("use_NN")),
    _torch_script_userobject(getUserObject<TorchScriptUserObject>("torch_script_userobject")),
    _input_tensor(torch::zeros(
        {1, 2},
        torch::TensorOptions().dtype(torch::kFloat64).device(_app.getLibtorchDevice()))),
    _property_prefix(getParam<MooseFunctorName>("property_prefix"))
{
    _properties.push_back(&declareGenericPropertyByName<Real, false>("ani_mu_t"));
    
    for (unsigned int i = 0; i < _mesh_dimension; ++i)
    {
      for (unsigned int j = 0; j < _mesh_dimension; ++j)
      {
        const auto name = _property_prefix + '_' + std::to_string(i) + std::to_string(j);
        _properties.push_back(&declareGenericPropertyByName<Real, false>(name));
      }
    }
}

void
TorchScriptTurbulentAnisotropyMaterial::initQpStatefulProperties()
{
  computeQpValues();
}

void
TorchScriptTurbulentAnisotropyMaterial::computeQpProperties()
{
  computeQpValues();
}

void 
TorchScriptTurbulentAnisotropyMaterial::arsm(const Real& eta1, const Real& eta2, Real& G1, Real& G2, Real& G3)
{
  /// SSG
  /*
  Real c_10 = 3.4;
  Real c_11 = 1.8;
  Real c_2  = 0.36;
  Real c_3  = 1.25;
  Real c_4 = 0.4;
  */

  /// LRR
  Real c_10 = 3.0;
  Real c_11 = 0.0;
  Real c_2  = 0.8;
  Real c_3  = 1.75;
  Real c_4 = 1.31;

  Real l_10 = c_10/2.0 - 1.0;
  Real l_11 = c_11 + 2.0;
  Real l_2 = c_2/2.0 - 2.0/3.0;
  Real l_3 = c_3/2.0 - 1.0;
  Real l_4 = c_4/2.0 - 1.0;

  Real p = - (2.0 * l_10) / (eta1 * l_11);
  Real q = (1.0 / std::pow(eta1 * l_11, 2.0)) * (std::pow(l_10, 2.0) + eta1 * l_11 * l_2 - (2.0/3.0) * eta1 * std::pow(l_3, 2.0) + 2.0 * eta2 * std::pow(l_4, 2.0));
  Real r = - (l_10 * l_2) / std::pow(eta1 * l_11, 2.0);

  Real a = (q - std::pow(p, 2.0) / 3.0);
  Real b = (1.0/27.0) * (2.0 * std::pow(p, 3.0) - 9.0 * p * q + 27.0 * r);
  Real d = std::pow(b, 2.0) / 4.0 + std::pow(a, 3.0) / 27.0;

  Real theta = std::acos(-(b / 2.0) / (std::sqrt(- std::pow(a, 3.0) / 27.0)));

  if (d > 0)
  {
    G1 = - (p / 3.0) + std::cbrt(-(b/2.0) + std::sqrt(d)) + std::cbrt(-(b/2.0) - std::sqrt(d));
  }
  else if (d <= 0 && b < 0)
  {
    G1 = - (p/3.0) + 2.0 * std::sqrt(-(a/3.0)) * std::cos(theta / 3.0);
  }
  else
  {
    G1 = - (p/3.0) + 2.0 * std::sqrt(-(a/3.0)) * std::cos(theta / 3.0 + (2.0 * libMesh::pi / 3.0));
  }

  G2 = - (l_4 * G1) / (l_10 - eta1 * l_11 * G1);
  G3 = (2.0 * l_3 * G1) / (l_10 - eta1 * l_11 * G1);
}


void
TorchScriptTurbulentAnisotropyMaterial::computeQpValues()
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

    const Real eta1 = std::min(std::max(sij.contract(sij), 1e-8), 1e8);
    const Real eta2 = std::min(std::max(rij.contract(rij), 1e-8), 1e8);

    double G1 = 0.0;
    double G2 = 0.0;
    double G3 = 0.0;

    if (_use_NN)
    {
    auto input_accessor = _input_tensor.accessor<Real, 2>();
    input_accessor[0][0] =  eta1;
    input_accessor[0][1] =  eta2;

    const auto output = _torch_script_userobject.evaluate(_input_tensor);
    const auto output_accessor = output.accessor<Real, 2>();

    G1 = output_accessor[0][0];
    G2 = output_accessor[0][1];
    G3 = output_accessor[0][2];
    }
    else
    {
        arsm(eta1, eta2, G1, G2, G3);
    }

    const Real rho = MetaPhysicL::raw_value(_rho(r,t));
    const TensorValue<Real> bij = rho * (G1 * sij * 0.0 + G2 * (sij * rij - rij * sij) + G3 * (sij * sij - 1./3. * I * sij.contract(sij)));
    const Real _ani_mu_t = - rho * k * G1 * timescale;

    bool irregular = false; 

    if (_debug || irregular)
    {
        _console << "-----------------------------------------" << std::endl;
        _console << "k: " << k << " eps: " << eps << std::endl;
        _console << "eta_1: " << eta1 << " eta_2: " << eta1 << std::endl;
        _console << "G1: " << G1 << " G2: " << G2 << " G3: " << G3 << std::endl;
    }
    if (irregular)
    {
        arsm(eta1, eta2, G1, G2, G3);
        _console << "ARSM G1: " << G1 << " G2: " << G2 << " G3: " << G3 << std::endl;
    }

    (*_properties[0])[_qp] = _ani_mu_t;
    for (unsigned int i = 0; i < _mesh_dimension; ++i)
    {
      for (unsigned int j = 0; j < _mesh_dimension; ++j)
      {
        const auto index = i * _mesh_dimension + j;
        (*_properties[index + 1])[_qp] = bij(i,j);
      }
    }
} 

#endif