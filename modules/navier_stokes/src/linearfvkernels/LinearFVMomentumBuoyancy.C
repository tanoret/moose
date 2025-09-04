//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "LinearFVMomentumBuoyancy.h"
#include "Assembly.h"
#include "SubProblem.h"
#include "NS.h"
#include "FEProblemBase.h"

registerMooseObject("NavierStokesApp", LinearFVMomentumBuoyancy);

InputParameters
LinearFVMomentumBuoyancy::validParams()
{
  InputParameters params = LinearFVElementalKernel::validParams();
  params.addClassDescription("Represents the buoyancy term in the Navier Stokes momentum "
                             "equations, added to the right hand side.");


  params.addRequiredParam<RealVectorValue>("gravity", "Gravitational acceleration vector.");
  params.addRequiredParam<MooseFunctorName>(NS::density, "The value for the density");
  params.addRequiredParam<Real>("reference_rho", "The value for the reference density");
  MooseEnum momentum_component("x=0 y=1 z=2");
  params.addRequiredParam<MooseEnum>(
      "momentum_component",
      momentum_component,
      "The component of the momentum equation that this kernel applies to.");

  return params;
}

LinearFVMomentumBuoyancy::LinearFVMomentumBuoyancy(const InputParameters & params)
  : LinearFVElementalKernel(params),
    _index(getParam<MooseEnum>("momentum_component")),
    _rho(getFunctor<Real>(NS::density)),
    _rho_0(getParam<Real>("reference_rho")),
    _gravity(getParam<RealVectorValue>("gravity"))
{
}

Real
LinearFVMomentumBuoyancy::computeMatrixContribution()
{
  return 0.0;
}

Real
LinearFVMomentumBuoyancy::computeRightHandSideContribution()
{
  const auto state = determineState();
  const auto & mesh = _subproblem.mesh();
  const auto coord_type       = mesh.getCoordSystem(_current_elem_info->elem()->subdomain_id());
  const auto rz_radial_coord  = mesh.getAxisymmetricRadialCoord();
  Real rhs = 0.0;
  auto action_functor = [this, state, &rhs](const Elem & /*elem*/,
                                            const Elem * /*neighbor*/,
                                            const FaceInfo * const fi,
                                            const Point & surface_vector,
                                            Real /*coord*/,
                                            const bool /*elem_has_info*/)
  {
    mooseAssert(fi, "FaceInfo is required.");
    // Use the FV functor to reconstruct ρ at the face with your chosen limiter.
    Moose::FaceArg face_arg = Moose::FaceArg{fi,
                                              Moose::FV::LimiterType::CentralDifference,
                                              true,
                                              /* correct_skewness */ false,
                                              this->_current_elem_info->elem(),
                                              nullptr};

    const Real rho_f       = _rho(face_arg, state);
    const Real rho_prime_f = rho_f - _rho_0;                       // Boussinesq/dynamic form
    rhs += rho_prime_f * _gravity(_index) * surface_vector(_index); // ∑_faces (\rho' g_i) (n_i A)_f
  };
  Moose::FV::loopOverElemFaceInfo(*_current_elem_info->elem(),
                                  mesh,
                                  action_functor,
                                  coord_type,
                                  rz_radial_coord);
  return rhs;
}
