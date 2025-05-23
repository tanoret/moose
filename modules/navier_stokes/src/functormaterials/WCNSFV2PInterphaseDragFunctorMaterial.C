//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "WCNSFV2PInterphaseDragFunctorMaterial.h"
#include "INSFVVelocityVariable.h"
#include "Function.h"
#include "NS.h"
#include "FVKernel.h"

registerMooseObject("NavierStokesApp", WCNSFV2PInterphaseDragFunctorMaterial);

InputParameters
WCNSFV2PInterphaseDragFunctorMaterial::validParams()
{
  InputParameters params = FunctorMaterial::validParams();
  params.addClassDescription("Computes the interphase drag coefficient for the two phase mixture model.");
  params.addRequiredCoupledVar("u", "The velocity in the x direction.");
  params.addCoupledVar("v", "The velocity in the y direction.");
  params.addCoupledVar("w", "The velocity in the z direction.");
  params.addRequiredCoupledVar("u_coupled", "The velocity in the x direction for the coupled phase.");
  params.addCoupledVar("v_coupled", "The velocity in the y direction for the coupled phase.");
  params.addCoupledVar("w_coupled", "The velocity in the z direction for the coupled phase.");

  params.addParam<MooseFunctorName>("fd_main", "Main phase void fraction - use only in problmes with >2 phases.");
  params.addRequiredParam<MooseFunctorName>("fd", "Dispersed phase fraction.");

  params.addRequiredParam<MooseFunctorName>(NS::density, "Phase density.");
  params.addRequiredParam<MooseFunctorName>(NS::density + "_coupled", "Coupled phase density.");

  params.addRequiredParam<MooseFunctorName>(NS::density, "Phase density.");
  params.addRequiredParam<MooseFunctorName>(NS::density + "_coupled", "Coupled phase density.");

  params.addRequiredParam<MooseFunctorName>(NS::mu, "Phase viscosity.");
  params.addRequiredParam<MooseFunctorName>(NS::mu + "_coupled", "Coupled phase viscosity.");

  params.addRequiredParam<MooseFunctorName>("rho_d", "Dispersed phase density.");

  params.addParam<MooseFunctorName>(
      "particle_diameter", 1.0, "Characteristic diameter of coupled phase, e.g., diameter of dispersed particles.");

  MooseEnum drag_formulation_type("mixture phase-specifc", "mixture");
  params.addParam<MooseEnum>(
      "drag_formulation_type", drag_formulation_type, "The method used for computing the Reynolds number in the drag function.");


  MooseEnum drag_model("schiller-naumann morsi-alexander", "schiller-nauman");
  params.addParam<MooseEnum>(
      "drag_model", drag_model, "The drag model used.");
  return params;
}

WCNSFV2PInterphaseDragFunctorMaterial::WCNSFV2PInterphaseDragFunctorMaterial(
    const InputParameters & params)
  : FunctorMaterial(params),
    _dim(_subproblem.mesh().dimension()),
    _u_var(params.isParamValid("u")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v", 0))
               : nullptr),
    _v_var(params.isParamValid("v")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v", 0))
               : nullptr),
    _w_var(params.isParamValid("w")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("w", 0))
               : nullptr),
    _u_var_coupled(params.isParamValid("u_coupled")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v", 0))
               : nullptr),
    _v_var_coupled(params.isParamValid("v_coupled")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("v_coupled", 0))
               : nullptr),
    _w_var_coupled(params.isParamValid("w_coupled")
               ? dynamic_cast<const INSFVVelocityVariable *>(getFieldVar("w_coupled", 0))
               : nullptr),
    _alpha_main(params.isParamValid("fd_main")
               ? dynamic_cast<const Moose::Functor<ADReal> *>(getFieldVar("fd_main", 0))
               : nullptr),
    _alpha(getFunctor<ADReal>("fd")),
    _rho(getFunctor<ADReal>(NS::density)),
    _rho_coupled(getFunctor<ADReal>(NS::density + "_coupled")),
    _mu(getFunctor<ADReal>(NS::mu)),
    _mu_coupled(getFunctor<ADReal>(NS::mu + "_coupled")),
    _particle_diameter(getFunctor<ADReal>("particle_diameter")),
    _drag_formulation_type(getParam<MooseEnum>("drag_formulation_type")),
    _drag_model(getParam<MooseEnum>("drag_model"))
{
  if (!_u_var)
    paramError("u", "The u velocity must be provided and be an INSFVVelocityVariable.");

  if (!_u_var_coupled)
    paramError("u_coupled", "The coupled phase u velocity must be provided and be an INSFVVelocityVariable.");

  if (_dim >= 2 && !_v_var)
    paramError("v",
               "In two or more dimensions, the v velocity must be supplied and it must be an "
               "INSFVVelocityVariable.");

  if (_dim >= 2 && !_v_var_coupled)
    paramError("v_coupled",
               "In two or more dimensions, the v_coupled velocity must be supplied and it must be an "
               "INSFVVelocityVariable.");

  if (_dim >= 3 && !_w_var)
    paramError("w",
               "In three-dimensions, the w velocity must be supplied and it must be an "
               "INSFVVelocityVariable.");

  if (_dim >= 3 && !_w_var)
    paramError("w",
               "In three-dimensions, the w velocity must be supplied and it must be an "
               "INSFVVelocityVariable.");

  addFunctorProperty<ADReal>(
    getParam<MooseFunctorName>("interphase_drag"),
    [this](const auto & r, const auto & t)
    {
        const auto alpha_2 = _alpha(r, t);
        const auto alpha_1 = _alpha_main ? (*_alpha_main)(r, t) : alpha_2;

        const auto relaxation_time = _rho_coupled(r, t) * Utility::pow<2>(_particle_diameter(r, t)) / (18.0 *_mu(r, t));

        ADRealVectorValue speed_diff = 0;
        speed_diff(0) = (*_u_var_coupled)(r, t) - (*_u_var)(r, t);
        if (_dim > 1)
        {
            speed_diff(1) = (*_v_var_coupled)(r, t) - (*_v_var)(r, t);
            if (_dim > 2)
                speed_diff(2) = (*_w_var_coupled)(r, t) - (*_w_var)(r, t);
        }
        const auto speed_dif_norm = speed_diff.norm();

        ADReal Re;
        if (_drag_formulation_type == "mixture")
        {
            const auto rho_mixture = alpha_1 * _rho(r, t) + alpha_2 * _rho_coupled(r, t);
            const auto mu_mixture = alpha_1 * _mu(r, t) + alpha_2 * _mu_coupled(r, t);
            Re = rho_mixture * speed_dif_norm * _particle_diameter(r, t) / mu_mixture;
        }
        else // (_drag_formulation_type == "phase-specific")
            Re = _rho(r, t) * speed_dif_norm * _particle_diameter(r, t) / _mu(r, t);

        ADReal CD;
        if (_drag_model == "schiller-naumann")
        {
            if (Re <= 1e3)
                CD = 24.0 * (1.0 + 0.15 * std::pow(Re, 0.687)) / (Re + 1e-20);
            else
                CD = 0.44;
        }
        else // _drag_model == "morsi-alexander"
        {
            ADReal a1, a2, a3;
            if (Re < 1e-1)
            {
                a1 = 0.0; a2 = 24.0; a3 = 0.0;
                Re = std::max(Re, 1e-20);
            }
            else if (Re >= 1e-1 && Re < 1e0)
            {
                a1 = 3.690; a2 = 22.730; a3 = 0.0903;
            }
            else if (Re >= 1e0 && Re < 1e1)
            {
                a1 = 1.222; a2 = 29.1667; a3 = 0.0903;
            }
            else if (Re >= 1e1 && Re < 1e2)
            {
                a1 = 0.6167; a2 = 46.50; a3 = -116.67;
            }
            else if (Re >= 1e2 && Re < 1e3)
            {
                a1 = 0.3644; a2 = 98.33; a3 = -2778.0;
            }
            else if (Re >= 1e3 && Re < 5e3)
            {
                a1 = 0.357; a2 = 148.62; a3 = -47.5e3;
            }
            else if (Re >= 5e3 && Re < 1e4)
            {
                a1 = 0.460; a2 = -490.546; a3 = 578.5e3;
            }
            else // (Re >= 1e4)
            {
                a1 = 0.5191; a2 = -1662.5; a3 = 5.4167e6;
            }

            CD = a1 + a2/Re + a3/Utility::pow<2>(Re);
        }

        const auto f = CD * Re / 24.0;

        return alpha_1 * alpha_2 * _rho_coupled(r, t) * f / relaxation_time;
    });
}
