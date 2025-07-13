//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "LinearFVFluxKernel.h"
#include "RhieChowMassFluxMultiPhase.h"
#include "LinearFVAdvectionDiffusionBC.h"

/**
 * An advection kernel that implements the advection term for the passive scalar transport equation.
 */
class LinearFVMultiPhaseFractionAdvection : public LinearFVFluxKernel
{
public:
  static InputParameters validParams();
  LinearFVMultiPhaseFractionAdvection(const InputParameters & params);

  virtual Real computeElemMatrixContribution() override;

  virtual Real computeNeighborMatrixContribution() override;

  virtual Real computeElemRightHandSideContribution() override;

  virtual Real computeNeighborRightHandSideContribution() override;

  virtual Real computeBoundaryMatrixContribution(const LinearFVBoundaryCondition & bc) override;

  virtual Real computeBoundaryRHSContribution(const LinearFVBoundaryCondition & bc) override;

  virtual void setupFaceData(const FaceInfo * face_info) override;

protected:
  /// Function to compute compression mass flux
  Real computeCompressionVelocityMassFlux(const Moose::FaceArg & face, const Real & grad_alpha);

  /// Function to compute the internal contribution from the compression velocity
  /// to the mass fluxes for the matrix term
  Real computeCompressionVelocityMassFluxMatrixContribution();

  /// Function to compute the internal contribution from the compression velocity
  /// to the mass fluxes for the RHS term
  Real computeCompressionVelocityMassFluxRHSContribution();

  /// The Rhie-Chow user object that provides us with the face velocity
  const RhieChowMassFluxMultiPhase & _mass_flux_provider;

  /// The dimension of the mesh
  const unsigned int _dim;

  /// The compression velocity control parameter
  const Real _c_alpha;

  /// The phase density used for compression velocity
  const Moose::Functor<Real> * _rho;

  /// x-velocity
  const Moose::Functor<Real> * _u_var;
  /// y-velocity
  const Moose::Functor<Real> * _v_var;
  /// z-velocity
  const Moose::Functor<Real> * _w_var;

  /// x-velocity-mixture
  const Moose::Functor<Real> * _u_var_mixture;
  /// y-velocity-mixture
  const Moose::Functor<Real> * _v_var_mixture;
  /// z-velocity-mixture
  const Moose::Functor<Real> * _w_var_mixture;

  /// Switch to enable/disable nonorthogonal correction in the stress term
  const bool _use_nonorthogonal_correction;

private:
  /// Container for the current advected interpolation coefficients on the face to make sure
  /// we don't compute it multiple times for different terms.
  std::pair<Real, Real> _advected_interp_coeffs;

  /// Container for the velocity on the face which will be reused in the advection term's
  /// matrix and right hand side contribution
  Real _total_adv_mass_face_flux;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;

  /// The limiter method
  Moose::FV::InterpMethod _limiter_method;

  /// Face argument for higher order face interpolation
  Moose::FaceArg _low_order_face;
  Moose::FaceArg _high_order_face;

  /// CMULES face limiter
  Real _lambda_f;
};
