//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVElectrophoresisSource.h"
#include "NS.h"
#include "NonlinearSystemBase.h"
#include "NavierStokesMethods.h"

registerMooseObject("ChemicalReactionsApp", FVElectrophoresisSource);

InputParameters
FVElectrophoresisSource::validParams()
{
  InputParameters params = FVFluxKernel::validParams();
  params.addClassDescription(
      "Computes the effective diffusion source term due to electropheresis.");

  params.addParam<Real>("F", 96485.3321, "Faraday Constant.");
  params.addParam<Real>("z", 1.0, "Number of electrons involved in the electrode reaction.");
  params.addParam<Real>("R", 8.31446261815324, "Universal gas constant.");
  params.addRequiredParam<MooseFunctorName>(NS::temperature, "The temperature field.");
  params.addRequiredParam<MooseFunctorName>("D", "The diffusion coefficient.");
  params.addRequiredParam<MooseFunctorName>("phi", "The electric potential.");

  params.set<unsigned short>("ghost_layers") = 2;

  return params;
}

FVElectrophoresisSource::FVElectrophoresisSource(const InputParameters & params)
  : FVFluxKernel(params),
    _dim(_subproblem.mesh().dimension()),
    _F(getParam<Real>("F")),
    _z(getParam<Real>("z")),
    _R(getParam<Real>("R")),
    _T(getFunctor<ADReal>(NS::temperature)),
    _D(getFunctor<ADReal>("D")),
    _phi(getFunctor<ADReal>("phi"))
{
}

ADReal
FVElectrophoresisSource::computeQpResidual()
{
  using namespace Moose::FV;
  const auto state = determineState();

  ADReal face_D;
  if (onBoundary(*_face_info))
    face_D = _D(makeCDFace(*_face_info), state);
  else
    Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                           face_D,
                           _D(elemArg(), state),
                           _D(neighborArg(), state),
                           *_face_info,
                           true);

  ADReal face_c;
  if (onBoundary(*_face_info))
    face_c = _var(makeCDFace(*_face_info), state);
  else
    Moose::FV::interpolate(Moose::FV::InterpMethod::Average,
                           face_c,
                           _var(elemArg(), state),
                           _var(neighborArg(), state),
                           *_face_info,
                           true);

  Moose::FaceArg face;
  const bool skewness_correction =
      (_var.faceInterpolationMethod() == Moose::FV::InterpMethod::SkewCorrectedAverage);
  if (onBoundary(*_face_info))
    face = singleSidedFaceArg();
  else
    face = makeCDFace(*_face_info, skewness_correction);

  const auto phi_grad = _phi.gradient(face, state) * _face_info->normal();

  return _F * _z / _R * (face_D * face_c * phi_grad);
}

void
FVElectrophoresisSource::computeResidual(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;

  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(std::make_pair(_var.number(), _var.sys().number()));
  auto r = MetaPhysicL::raw_value(fi.faceArea() * fi.faceCoord() * computeQpResidual());

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    // residual contribution of this kernel to the elem element
    prepareVectorTag(_assembly, _var.number());
    _local_re(0) = r / MetaPhysicL::raw_value(_T(elemArg(), determineState()));
    accumulateTaggedLocalResidual();
  }
  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    // residual contribution of this kernel to the neighbor element
    prepareVectorTagNeighbor(_assembly, _var.number());
    _local_re(0) = -r / MetaPhysicL::raw_value(_T(neighborArg(), determineState()));
    accumulateTaggedLocalResidual();
  }
}

void
FVElectrophoresisSource::computeJacobian(const FaceInfo & fi)
{
  if (skipForBoundary(fi))
    return;

  _face_info = &fi;
  _normal = fi.normal();
  _face_type = fi.faceType(std::make_pair(_var.number(), _var.sys().number()));
  const ADReal r = fi.faceArea() * fi.faceCoord() * computeQpResidual();

  if (_face_type == FaceInfo::VarFaceNeighbors::ELEM ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    mooseAssert(_var.dofIndices().size() == 1, "We're currently built to use CONSTANT MONOMIALS");

    addResidualsAndJacobian(_assembly,
                            std::array<ADReal, 1>{{r / _T(elemArg(), determineState())}},
                            _var.dofIndices(),
                            _var.scalingFactor());
  }

  if (_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR ||
      _face_type == FaceInfo::VarFaceNeighbors::BOTH)
  {
    mooseAssert((_face_type == FaceInfo::VarFaceNeighbors::NEIGHBOR) ==
                    (_var.dofIndices().size() == 0),
                "If the variable is only defined on the neighbor hand side of the face, then that "
                "means it should have no dof indices on the elem element. Conversely if "
                "the variable is defined on both sides of the face, then it should have a non-zero "
                "number of degrees of freedom on the elem element");

    // We switch the sign for the neighbor residual
    ADReal neighbor_r = -r / _T(neighborArg(), determineState());

    mooseAssert(_var.dofIndicesNeighbor().size() == 1,
                "We're currently built to use CONSTANT MONOMIALS");

    addResidualsAndJacobian(_assembly,
                            std::array<ADReal, 1>{{neighbor_r}},
                            _var.dofIndicesNeighbor(),
                            _var.scalingFactor());
  }
}
