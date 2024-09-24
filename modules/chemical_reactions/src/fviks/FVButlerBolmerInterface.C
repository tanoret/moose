//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "FVButlerBolmerInterface.h"
#include "NS.h"

registerMooseObject("ChemicalReactionsApp", FVButlerBolmerInterface);

InputParameters
FVButlerBolmerInterface::validParams()
{
  InputParameters params = FVInterfaceKernel::validParams();
  params.addClassDescription("Computes the residual for coupled redox reactions across interfaces "
                             "using the Butler-Volmer formulation.");

  // Background parameters
  params.addParam<Real>("F", 96485.3321, "Faraday Constant.");
  params.addParam<Real>("z", 1.0, "Number of electrons involved in the electrode reaction.");
  params.addParam<Real>("R", 8.31446261815324, "Universal gas constant.");
  params.addRequiredParam<MooseFunctorName>(NS::temperature, "The temperature field.");
  params.addRequiredParam<MooseFunctorName>("phi", "The electric potential.");

  // Tracked Species Parameters
  params.addRequiredParam<MooseFunctorName>("c",
                                            "The main concentration field in the fluid domain.");
  params.addParam<MooseFunctorName>("c_solid", "The main concentration field in the solid domain.");
  params.addRequiredParam<Real>("c_E", "Electrode potential for the main tracked species.");
  params.addParam<Real>(
      "c_dE_dT",
      0.0,
      "Temperature derivative of electrode potential for the main tracked species.");
  params.addParam<Real>("T_ref", 300.0, "Reference electrode temperature.");
  params.addRequiredParam<unsigned int>(
      "c_Z", "Electrode transfered charge for the main tracked species.");
  params.addParam<Real>("c_G", 0.0, "Gibbs formation energy for the main tracked species.");

  // Formulation used for Butler-Volmer current
  MooseEnum interface_formulation_enum("standard multi-species", "standard");
  params.addParam<MooseEnum>("interface_formulation",
                             interface_formulation_enum,
                             "Formulation to use for the Butler-Volmer current.");

  /// Parameters spcific to standard half-cell formulation
  params.addParam<Real>("c_k0", "External reaction constant.");
  params.addParam<Real>("c_alpha", 0.5, "Charge transfer coefficient for half-cell reaction.");

  // Coupled Species Parameters for Multi-Species Formulation
  params.addParam<std::vector<MooseFunctorName>>(
      "coupled_c",
      {},
      "The corresponding names of the functors for the couples "
      "species for electro-chemistry modeling.");
  params.addParam<std::vector<MooseFunctorName>>(
      "coupled_c_solid",
      "The corresponding names of the functors for the couples "
      "species for electro-chemistry modeling in the solid domain.");
  params.addParam<std::vector<Real>>("coupled_E", {0.0}, "Electrode potential of coupled species.");
  params.addParam<std::vector<Real>>(
      "coupled_dE_dT", "Temperature derivative of electrode potential of coupled species.");
  params.addParam<std::vector<unsigned int>>(
      "coupled_Z", {1}, "Electrode transfered charge of coupled species.");
  params.addParam<std::vector<Real>>("coupled_G", "Gibbs formation energy of coupled species.");
  params.addParam<std::vector<Real>>("coupled_k0", "Reaction constants of coupled species.");
  params.addParam<std::vector<Real>>(
      "coupled_alpha", "Anodic charge transfer coefficient of the reaction with couoked species.");

  // Interface configuration transfers
  params.addParam<Real>(
      "bulk_distance", -1, "The distance to the bulk for evaluating the fluid bulk temperature");
  params.addParam<bool>("wall_cell_is_bulk",
                        false,
                        "Use the wall cell centroid temperature for the fluid bulk temperature");
  return params;
}

FVButlerBolmerInterface::FVButlerBolmerInterface(const InputParameters & params)
  : FVInterfaceKernel(params),
    _F(getParam<Real>("F")),
    _z(getParam<Real>("z")),
    _R(getParam<Real>("R")),
    _T(getFunctor<ADReal>(NS::temperature)),
    _phi(getFunctor<ADReal>("phi")),
    _c(getFunctor<ADReal>("c")),
    _c_solid(params.isParamValid("c_solid") ? &getFunctor<ADReal>("c_solid")
                                            : &getFunctor<ADReal>("c")),
    _c_E(getParam<Real>("c_E")),
    _c_dE_dT(getParam<Real>("c_dE_dT")),
    _T_ref(getParam<Real>("make ")),
    _c_Z(getParam<unsigned int>("c_Z")),
    _c_G(getParam<Real>("c_G")),
    _interface_formulation(getParam<MooseEnum>("interface_formulation")),
    _c_k0(params.isParamValid("c_k0") ? &getParam<Real>("c_k0") : nullptr),
    _c_alpha(getParam<Real>("c_alpha")),
    _coupled_c(getParam<std::vector<MooseFunctorName>>("coupled_c")),
    _coupled_c_solid(params.isParamValid("coupled_c_solid")
                         ? &getParam<std::vector<MooseFunctorName>>("coupled_c_solid")
                         : &getParam<std::vector<MooseFunctorName>>("coupled_c")),
    _coupled_E(getParam<std::vector<Real>>("coupled_E")),
    _coupled_dE_dT(params.isParamValid("coupled_dE_dT")
                       ? &getParam<std::vector<Real>>("coupled_dE_dT")
                       : nullptr),
    _coupled_Z(getParam<std::vector<unsigned int>>("coupled_Z")),
    _coupled_G(params.isParamValid("coupled_G") ? &getParam<std::vector<Real>>("coupled_G")
                                                : nullptr),
    _coupled_k0(params.isParamValid("coupled_k0") ? &getParam<std::vector<Real>>("coupled_k0")
                                                  : nullptr),
    _coupled_alpha(params.isParamValid("coupled_alpha")
                       ? &getParam<std::vector<Real>>("coupled_alpha")
                       : nullptr),
    _bulk_distance(getParam<Real>("bulk_distance")),
    _use_wall_cell(getParam<bool>("wall_cell_is_bulk")),
    _pl(mesh().getPointLocator()),
    _var1_is_fluid("wraps_" + var1().name() == _c.functorName())
{
  if (!_use_wall_cell && (_bulk_distance < 0))
    mooseError(
        "The bulk distance should be specified or 'wall_cell_is_bulk' should be set to true for "
        "the FVTwoVarConvectionCorrelationInterface");

  if (_interface_formulation == "multi-species")
  {

    // Getting the functors for coupled species
    _coupled_c_functors.resize(_coupled_c.size());
    for (const auto i : make_range(_coupled_c.size()))
      _coupled_c_functors[i] = &getFunctor<GenericReal<true>>((*_coupled_c_solid)[i]);

    if ((*_coupled_c_solid).size() != _coupled_c.size())
      paramError("coupled_c_solid",
                 "The number of spceified coupled species in the solid domain is different than "
                 "the number of "
                 "coupled species 'c'");

    _coupled_c_solid_functors.resize((*_coupled_c_solid).size());
    for (const auto i : make_range((*_coupled_c_solid).size()))
      _coupled_c_solid_functors[i] = &getFunctor<GenericReal<true>>((*_coupled_c_solid)[i]);

    if (_coupled_E.size() != _coupled_c.size())
      paramError("coupled_E",
                 "The number of spceified electrode reactions is different than the number of "
                 "coupled species 'c'");

    if (_coupled_Z.size() != _coupled_c.size())
      paramError("coupled_Z",
                 "The number of spceified transfered charges in electrode reactions is different "
                 "than the number of "
                 "coupled species 'c'");

    if (_coupled_dE_dT)
    {
      if (_coupled_dE_dT->size() != _coupled_c.size())
        paramError("coupled_dE_dT",
                   "The number of specified temperature derivatives of electrode potentials is "
                   "different than the number of "
                   "coupled species 'c'");
    }
    else
    {
      _coupled_dE_dT = new std::vector<Real>(_coupled_c.size(), 0.0);
    }

    if (_coupled_G)
    {
      if (_coupled_G->size() != _coupled_c.size())
        paramError(
            "coupled_G",
            "The number of specified Gibbs formation energies is different than the number of "
            "coupled species 'c'");
    }
    else
    {
      _coupled_G = new std::vector<Real>(_coupled_c.size(), 0.0);
    }

    if (_coupled_k0)
      if (_coupled_k0->size() != _coupled_c.size())
        paramError(
            "coupled_k0",
            "The number of specified Gibbs formation energies is different than the number of "
            "coupled species 'c'");

    if (_coupled_alpha)
    {
      if (_coupled_alpha->size() != _coupled_c.size())
        paramError(
            "coupled_alpha",
            "The number of specified Gibbs formation energies is different than the number of "
            "coupled species 'c'");
    }
    else
    {
      _coupled_alpha = new std::vector<Real>(_coupled_c.size(), 0.5);
    }
  }
}

ADReal
FVButlerBolmerInterface::computeQpResidual()
{
  // If variable1 is fluid and variable 1 is on elem or
  // if variable2 is fluid and variable 2 is on elem
  // the fluid element will be elem otherwise it is the neighbor
  const Elem * elem_on_fluid_side =
      (elemIsOne() && _var1_is_fluid) || (!elemIsOne() && !_var1_is_fluid)
          ? &_face_info->elem()
          : _face_info->neighborPtr();

  const Elem * bulk_elem;
  const auto state = determineState();
  if (!_use_wall_cell)
  {
    Point p = _face_info->faceCentroid();
    Point du = Point(MetaPhysicL::raw_value(_normal));
    du *= _bulk_distance;
    // The normal always points outwards from the elem (towards the neighbor)
    if (elem_on_fluid_side == &_face_info->elem())
      p -= du;
    else
      p += du;
    bulk_elem = (*_pl)(p);
  }
  else
    bulk_elem = elem_on_fluid_side;

  mooseAssert(bulk_elem,
              "The element at bulk_distance from the wall was not found in the mesh. "
              "Increase the number of ghost layers with the 'ghost_layers' parameter.");
  mooseAssert((_var1_is_fluid ? var1() : var2()).hasBlocks(bulk_elem->subdomain_id()),
              "The fluid temperature is not defined at bulk_distance from the wall.");

  const auto fluid_side = singleSidedFaceArg(_var1_is_fluid ? var1() : var2(), _face_info);
  const auto solid_side = singleSidedFaceArg(_var1_is_fluid ? var2() : var1(), _face_info);

  const auto bulk_elem_arg = makeElemArg(bulk_elem);

  // Build up Butler-Volmer current
  ADReal c_current = 0.0;
  const unsigned int loop_length =
      (_interface_formulation == "multi-species") ? _coupled_c.size() : 1;
  for (unsigned int i = 0; i < loop_length; ++i)
  {
    const auto Z = (_c_Z / std::gcd(_c_Z, _coupled_Z[i])) * _coupled_Z[i];
    const auto T_bulk = _T(bulk_elem_arg, state);
    const auto ZFoRT = Z * _F / (_R * T_bulk);
    const auto background_potential = _phi(bulk_elem_arg, state);

    const Real loc_alpha = (_interface_formulation == "standard") ? _c_alpha : (*_coupled_alpha)[i];
    const Real coupled_E = (_interface_formulation == "standard") ? 0.0 : _coupled_E[i];
    const Real coupled_dE_dT = (_interface_formulation == "standard") ? 0.0 : (*_coupled_dE_dT)[i];

    const auto electrode_potential =
        (coupled_E + coupled_dE_dT * (T_bulk - _T_ref) - _c_E -
         _c_dE_dT * (T_bulk - _T_ref)); //+ std::log(c_int / c_coupled_int) / ZFoRT;
    const auto effective_potential = electrode_potential - background_potential;

    ADReal oxidizing_partial_current, reducing_partial_current;

    if (_interface_formulation == "standard")
    {
      const auto c_int_fluid = _c(fluid_side, state);
      const auto c_int_solid = (*_c_solid)(solid_side, state);
      oxidizing_partial_current =
          c_int_solid * std::exp(ZFoRT * (1.0 - loc_alpha) * effective_potential);
      reducing_partial_current = c_int_fluid * std::exp(-ZFoRT * loc_alpha * effective_potential);
    }
    else
    {
      const auto c_int_fluid_tracked = std::max(_c(fluid_side, state), 1e-10);
      const auto c_int_solid_tracked = std::max((*_c_solid)(solid_side, state), 1e-10);
      const auto species_oxidizing_partial_current =
          c_int_solid_tracked * std::exp(ZFoRT * (1.0 - loc_alpha) * effective_potential);
      const auto species_reducing_partial_current =
          c_int_fluid_tracked * std::exp(-ZFoRT * loc_alpha * effective_potential);

      const auto c_int_fluid_reac = std::max((*_coupled_c_functors[i])(fluid_side, state), 1e-10);
      const auto c_int_solid_reac =
          std::max((*_coupled_c_solid_functors[i])(solid_side, state), 1e-10);
      const auto coupled_oxidizing_partial_current =
          c_int_solid_reac * std::exp(ZFoRT * (1.0 - loc_alpha) * effective_potential);
      const auto coupled_reducing_partial_current =
          c_int_fluid_reac * std::exp(-ZFoRT * loc_alpha * effective_potential);

      oxidizing_partial_current =
          species_oxidizing_partial_current * 0.0 + coupled_oxidizing_partial_current;
      reducing_partial_current =
          species_reducing_partial_current * 0.0 + coupled_reducing_partial_current;
    }

    ADReal k0;
    if (_interface_formulation == "standard")
    {
      if (_c_k0)
        k0 = (*_c_k0) * std::pow(std::max((*_c_solid)(solid_side, state), 1e-10), (1 - _c_alpha)) *
             std::pow(std::max(_c(fluid_side, state), 1e-10), _c_alpha);
      else
        k0 = std::exp(-_c_G / (_R * T_bulk));
    }
    else
    {
      if (_coupled_k0)
        k0 = (*_coupled_k0)[i] *
             std::pow(std::max((*_c_solid)(solid_side, state), 1e-10), (1 - _c_alpha)) *
             std::pow(std::max(_c(fluid_side, state), 1e-10), _c_alpha);
      else
        k0 = std::exp(-(_c_G - (*_coupled_G)[i]) / (_R * T_bulk));
    }

    c_current += k0 * (oxidizing_partial_current - reducing_partial_current);
  }

  auto multipler =
      _normal * (_face_info->faceCentroid() - bulk_elem->vertex_average()) > 0 ? 1 : -1;

  return -multipler * c_current;
}
