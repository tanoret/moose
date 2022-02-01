#pragma once

#include "INSFVFluxKernel.h"
#include "INSFVMomentumResidualObject.h"

/// INSFVMomentumRegularization implements a standard diffusion term:
///
///     - strong form: \nabla \cdot k \nabla u
///
///     - weak form: \int_{A} k \nabla u \cdot \vec{n} dA
///
/// It uses/requests a material property named "_filter_scaling" for k. An average of
/// the elem and neighbor k-values (which should be face-values) is used to
/// compute k on the face. Cross-diffusion correction factors are currently not
/// implemented for the "grad_u*n" term.
class INSFVMomentumRegularization : public INSFVFluxKernel
{
public:
  static InputParameters validParams();
  INSFVMomentumRegularization(const InputParameters & params);

  // This object neither contributes to the A coefficients nor to the B (source) coefficients
  using INSFVFluxKernel::gatherRCData;
  void gatherRCData(const FaceInfo &) override {}

protected:
  virtual ADReal computeQpResidual() override;

  const Moose::Functor<ADReal> & _filter_scaling;
};
