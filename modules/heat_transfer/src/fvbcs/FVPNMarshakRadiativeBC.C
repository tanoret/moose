//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVPNMarshakRadiativeBC.h"
#include "DenseMatrix.h"

registerMooseObject("HeatTransferApp", FVPNMarshakRadiativeBC);

InputParameters
FVPNMarshakRadiativeBC::validParams()
{
  InputParameters params = FVDirichletBCBase::validParams();
  params.addClassDescription("Adds the Marshak boundary conditions for a radative flux with order N");

    params.addRequiredParam<unsigned int>("n", "N-th Order of the Equation.");

    params.addParam<MooseFunctorName>("phi_0", "Order 0 of the radiative flux.");
    params.addParam<MooseFunctorName>("phi_2", "Order 2 of the radiative flux.");
    params.addParam<MooseFunctorName>("phi_4", "Order 4 of the radiative flux.");
    params.addParam<MooseFunctorName>("phi_6", "Order 6 of the radiative flux.");

    params.addParam<MooseFunctorName>("boundary_source", "Diffuse source at the boundary.");

  return params;
}

FVPNMarshakRadiativeBC::FVPNMarshakRadiativeBC(const InputParameters & params)
  : FVDirichletBCBase(params),
    _n(getParam<unsigned int>("n")),
    _phi_0(params.isParamValid("phi_0") ? &(getFunctor<ADReal>("phi_0")) : nullptr),
    _phi_2(params.isParamValid("phi_2") ? &(getFunctor<ADReal>("phi_2")) : nullptr),
    _phi_4(params.isParamValid("phi_4") ? &(getFunctor<ADReal>("phi_4")) : nullptr),
    _phi_6(params.isParamValid("phi_6") ? &(getFunctor<ADReal>("phi_6")) : nullptr),
    _boundary_source(params.isParamValid("boundary_source") ? &(getFunctor<ADReal>("boundary_source")) : nullptr)
{
    _marshak_inlet_matrix = new DenseMatrix<Real>(4,4);

    (*_marshak_inlet_matrix)(0, 0) = 1./2.;
    (*_marshak_inlet_matrix)(0, 1) = -1./8.;
    (*_marshak_inlet_matrix)(0, 2) = 1./16.;
    (*_marshak_inlet_matrix)(0, 3) = -5./128.;

    (*_marshak_inlet_matrix)(1, 0) = -1./8.;
    (*_marshak_inlet_matrix)(1, 1) = 7./24.;
    (*_marshak_inlet_matrix)(1, 2) = -41./384.;
    (*_marshak_inlet_matrix)(1, 3) = 1./16.;

    (*_marshak_inlet_matrix)(2, 0) = 1./16.;
    (*_marshak_inlet_matrix)(2, 1) = -41./384.;
    (*_marshak_inlet_matrix)(2, 2) = 407./1920.;
    (*_marshak_inlet_matrix)(2, 3) = -233./2560.;

    (*_marshak_inlet_matrix)(3, 0) = -5./128.;
    (*_marshak_inlet_matrix)(3, 1) = 1./16.;
    (*_marshak_inlet_matrix)(3, 2) = -233./2560.;
    (*_marshak_inlet_matrix)(3, 3) = 3023./17920.;

    _N = 0;
    if (_phi_2)
        _N += 1;
    if (_phi_4)
        _N += 1;
    if (_phi_6)
        _N += 1;  

    _phi_vector.push_back(_phi_0);
    _phi_vector.push_back(_phi_2);
    _phi_vector.push_back(_phi_4);
    _phi_vector.push_back(_phi_6);

}

ADReal
FVPNMarshakRadiativeBC::boundaryValue(const FaceInfo & fi, const Moose::StateArg & state) const
{
    const auto boundary_face = singleSidedFaceArg(&fi);

    const unsigned int n_matrix = _n / 2;

    const auto diag_coef = (*_marshak_inlet_matrix)(n_matrix, n_matrix);

    ADReal rhs(0.0);
    for (unsigned int i = 0; i < _N; ++i)
    {
        if (i != n_matrix)
            rhs += (*_marshak_inlet_matrix)(n_matrix, i) * (*_phi_vector[i])(boundary_face, state);
    }

    if (_boundary_source && _n == 0)
        rhs += (*_boundary_source)(boundary_face, state);
  
    return rhs / diag_coef;
}
