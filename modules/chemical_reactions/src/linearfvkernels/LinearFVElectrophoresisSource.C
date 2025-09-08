//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVElectrophoresisSource.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "FEProblemBase.h"

registerMooseObject("ChemicalReactionsApp", LinearFVElectrophoresisSource);

InputParameters
LinearFVElectrophoresisSource::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();

  params.addClassDescription(
      "Computes the effective diffusion source term due to electropheresis.");

  params.addParam<Real>("F", 96485.3321, "Faraday Constant.");
  params.addParam<Real>("z", 1.0, "Number of electrons involved in the electrode reaction.");
  params.addParam<Real>("R", 8.31446261815324, "Universal gas constant.");
  params.addRequiredParam<MooseFunctorName>("T", "The temperature field.");
  params.addRequiredParam<MooseFunctorName>("D", "The diffusion coefficient.");
  params.addParam<VariableName>("phi",
                                "The electric potential variable whose gradient should be used.");
  
  return params;
}

LinearFVElectrophoresisSource::LinearFVElectrophoresisSource(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _F(getParam<Real>("F")),
    _z(getParam<Real>("z")),
    _R(getParam<Real>("R")),
    _T(getFunctor<Real>("T")),
    _D(getFunctor<Real>("D")),
    _electric_potential_var(getelectric_potentialVariable("phi")),
    _electric_potential_gradient(_electric_potential_var.sys().gradientContainer()),
    _electric_potential_var_num(_electric_potential_var.number()),
    _electric_potential_sys_num(_electric_potential_var.sys().number())
{
  _electric_potential_var.computeCellGradients();
}

MooseLinearVariableFV<Real> &
LinearFVElectrophoresisSource::getelectric_potentialVariable(const std::string & vname)
{
  auto * ptr = dynamic_cast<MooseLinearVariableFV<Real> *>(
      &_fe_problem.getVariable(_tid, getParam<VariableName>(vname)));

  if (!ptr)
    paramError("phi", "The electric_potential variable should be of type MooseLinearVariableFVReal!");

  return *ptr;
}

Real
LinearFVElectrophoresisSource::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVElectrophoresisSource::computeRightHandSideContribution()
{
  const auto elem_arg = makeElemArg(_current_elem_info->elem());
  const auto coord_type = _subproblem.mesh().getCoordSystem(elem_arg.elem->subdomain_id());
  const auto rz_radial_coord = _subproblem.mesh().getAxisymmetricRadialCoord();

  Real grad_term = 0.0;
  RealVectorValue gradient_elem(0.0);

  std::vector<RealVectorValue> A_sys;
  std::vector<Real> rhs_sys;

  auto action_functor = [this,
                         &elem_arg,
                         &grad_term,
                         &gradient_elem]
                                (const Elem & /*elem*/,
                                const Elem * /*neighbor*/,
                                const FaceInfo * const fi,
                                const Point & surface_vector,
                                Real /*coord*/,
                                const bool /*elem_has_info*/)
  {
      mooseAssert(fi, "We need a FaceInfo for this action_functor");

      Moose::FaceArg face_arg = Moose::FaceArg{fi,
                                               Moose::FV::LimiterType::CentralDifference,
                                               true,
                                               /* correct_skewness */ false,
                                               this->_current_elem_info->elem(),
                                               nullptr};

      const auto grad_electric_potential = MetaPhysicL::raw_value(this->_electric_potential_var.gradient(face_arg, this->determineState()));
      constexpr Real tiny = 1.0e-14;
      auto grad_term_contrib = grad_electric_potential * surface_vector;
      grad_term_contrib *= _D(face_arg, this->determineState()) * 
                          Utility::pow<2>(_z * _F) /
                          std::max(_R * _T(face_arg, this->determineState()), tiny);
      grad_term += grad_term_contrib;
  };

  Moose::FV::loopOverElemFaceInfo(
      *_current_elem_info->elem(), _subproblem.mesh(), action_functor, coord_type, rz_radial_coord);

  return grad_term;
}
