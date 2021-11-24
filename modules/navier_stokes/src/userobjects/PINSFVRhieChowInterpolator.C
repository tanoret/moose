//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "PINSFVRhieChowInterpolator.h"
#include "Reconstructions.h"
#include "NS.h"

registerMooseObject("NavierStokesApp", PINSFVRhieChowInterpolator);

InputParameters
PINSFVRhieChowInterpolator::validParams()
{
  auto params = INSFVRhieChowInterpolator::validParams();
  params.addClassDescription("Performs interpolations and reconstructions of body forces and "
                             "porosity and computes the Rhie-Chow face velocities.");
  params.addRequiredParam<MooseFunctorName>(NS::porosity, "The porosity");
  params.addParam<unsigned short>(
      "reconstructions", 0, "The number of reconstructions to perform on the porosity");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC,
      [](const InputParameters & obj_params, InputParameters & rm_params) {
        rm_params.set<unsigned short>("layers") = obj_params.get<unsigned short>("reconstructions");
        rm_params.set<bool>("use_displaced_mesh") = obj_params.get<bool>("use_displaced_mesh");
      });
  params.addParam<bool>(
      "smooth_porosity", false, "Whether the porosity field is smooth or has discontinuities");
  return params;
}

PINSFVRhieChowInterpolator::PINSFVRhieChowInterpolator(const InputParameters & params)
  : INSFVRhieChowInterpolator(params),
    _eps(const_cast<Moose::Functor<ADReal> &>(getFunctor<ADReal>(NS::porosity))),
    _epss(libMesh::n_threads(), nullptr),
    _rec(getParam<unsigned short>("reconstructions")),
    _reconstructed_eps(_moose_mesh, true),
    _smooth_porosity(getParam<bool>("smooth_porosity"))
{
  if (_rec && _eps.wrapsType<MooseVariableBase>())
    paramError(
        NS::porosity,
        "If we are reconstructing porosity, then the input porosity to this user object cannot "
        "be a Moose variable. There are issues with reconstructing Moose variables: 1) initial "
        "conditions are run after use object initial setup 2) reconstructing from a variable "
        "requires ghosting the solution vectors 3) it's difficult to restrict the face "
        "informations we evaluate reconstructions on such that we never query an algebraically "
        "remote element due to things like two-term extrapolated boundary faces which trigger "
        "gradient evaluations which trigger neighbor element evaluation");

  const auto porosity_name = deduceFunctorName(NS::porosity);

  for (const auto tid : make_range(libMesh::n_threads()))
    _epss[tid] = &UserObject::_subproblem.getFunctor<ADReal>(porosity_name, tid, name());
}

void
PINSFVRhieChowInterpolator::meshChanged()
{
  INSFVRhieChowInterpolator::meshChanged();
  pinsfvSetup();
}

void
PINSFVRhieChowInterpolator::residualSetup()
{
  INSFVRhieChowInterpolator::residualSetup();
  if (!_initial_setup_done)
    pinsfvSetup();

  _initial_setup_done = true;
}

void
PINSFVRhieChowInterpolator::pinsfvSetup()
{
  if (!_rec)
    return;

  const auto & all_fi = _moose_mesh.allFaceInfo();
  _geometric_fi.reserve(all_fi.size());

  for (const auto & fi : all_fi)
    if (isFaceGeometricallyRelevant(fi))
      _geometric_fi.push_back(&fi);

  _geometric_fi.shrink_to_fit();

  const auto saved_do_derivatives = ADReal::do_derivatives;
  ADReal::do_derivatives = true;
  _reconstructed_eps.mapFilled(false);

  Moose::FV::reconstruct(_reconstructed_eps, _eps, _rec, false, false, _geometric_fi, *this);
  ADReal::do_derivatives = saved_do_derivatives;
  _reconstructed_eps.mapFilled(true);

  _eps.assign(_reconstructed_eps);

  const auto porosity_name = deduceFunctorName(NS::porosity);
  for (const auto tid : make_range((unsigned int)(1), libMesh::n_threads()))
  {
    auto & other_epss = const_cast<Moose::Functor<ADReal> &>(
        UserObject::_subproblem.getFunctor<ADReal>(porosity_name, tid, name()));
    other_epss.assign(_reconstructed_eps);
  }
}

VectorValue<ADReal>
PINSFVRhieChowInterpolator::getVelocity(const Moose::FV::InterpMethod m,
                                        const FaceInfo & fi,
                                        const THREAD_ID tid) const
{
  const Elem * const elem = &fi.elem();
  const Elem * const neighbor = fi.neighborPtr();
  auto & vel = *_vel[tid];
  auto & eps = *_epss[tid];
  auto & p = *_ps[tid];

  if (Moose::FV::onBoundary(*this, fi))
  {
    const auto sub_id =
        hasBlocks(elem->subdomain_id()) ? elem->subdomain_id() : neighbor->subdomain_id();
    const auto boundary_face =
        std::make_tuple(&fi, Moose::FV::LimiterType::CentralDifference, true, sub_id);
    return vel(boundary_face);
  }

  VectorValue<ADReal> velocity;

  const auto elem_face = Moose::FV::elemFromFace(*this, fi);
  const auto neighbor_face = Moose::FV::neighborFromFace(*this, fi);

  Moose::FV::interpolate(
      Moose::FV::InterpMethod::Average, velocity, vel(elem_face), vel(neighbor_face), fi, true);

  if (m == Moose::FV::InterpMethod::Average)
    return velocity;

  // Avoid computing a pressure gradient near porosity jumps. We may remove this soon; see
  // https://github.com/idaholab/moose/issues/19425
  if (!_smooth_porosity)
    if (MetaPhysicL::raw_value(eps.gradient(elem)).norm() > 1e-12 ||
        MetaPhysicL::raw_value(eps.gradient(neighbor)).norm() > 1e-12)
      return velocity;

  mooseAssert(neighbor && this->hasBlocks(neighbor->subdomain_id()),
              "We should be on an internal face...");

  // Get pressure gradient. This is the uncorrected gradient plus a correction from cell centroid
  // values on either side of the face
  const VectorValue<ADReal> & grad_p = p.adGradSln(fi);

  // Get uncorrected pressure gradient. This will use the element centroid gradient if we are
  // along a boundary face
  const VectorValue<ADReal> & unc_grad_p = p.uncorrectedAdGradSln(fi);

  const Point & elem_centroid = fi.elemCentroid();
  const Point & neighbor_centroid = fi.neighborCentroid();
  Real elem_volume = fi.elemVolume();
  Real neighbor_volume = fi.neighborVolume();

  // Now we need to perform the computations of D
  const VectorValue<ADReal> & elem_a = libmesh_map_find(_a, elem->id());

  mooseAssert(UserObject::_subproblem.getCoordSystem(elem->subdomain_id()) ==
                  UserObject::_subproblem.getCoordSystem(neighbor->subdomain_id()),
              "Coordinate systems must be the same between the two elements");

  Real coord;
  coordTransformFactor(UserObject::_subproblem, elem->subdomain_id(), elem_centroid, coord);

  elem_volume *= coord;

  VectorValue<ADReal> elem_D = 0;
  for (const auto i : make_range(_dim))
  {
    mooseAssert(elem_a(i).value() != 0, "We should not be dividing by zero");
    elem_D(i) = elem_volume / elem_a(i);
  }

  VectorValue<ADReal> face_D;

  const VectorValue<ADReal> & neighbor_a = libmesh_map_find(_a, neighbor->id());

  coordTransformFactor(UserObject::_subproblem, neighbor->subdomain_id(), neighbor_centroid, coord);
  neighbor_volume *= coord;

  VectorValue<ADReal> neighbor_D = 0;
  for (const auto i : make_range(_dim))
  {
    mooseAssert(neighbor_a(i).value() != 0, "We should not be dividing by zero");
    neighbor_D(i) = neighbor_volume / neighbor_a(i);
  }
  Moose::FV::interpolate(Moose::FV::InterpMethod::Average, face_D, elem_D, neighbor_D, fi, true);

  // evaluate face porosity, see (18) in Hanimann 2021 or (11) in Nordlund 2016
  const auto face_eps = eps(Moose::FV::makeCDFace(fi, Moose::FV::faceArgSubdomains(*this, fi)));

  const auto b1 = _b(fi);
  const auto b3 = _b2(fi);

  // perform the pressure correction
  for (const auto i : make_range(_dim))
    velocity(i) +=
        -face_D(i) * face_eps * (grad_p(i) - unc_grad_p(i)) + face_D(i) * (b1(i) - b3(i));

  return velocity;
}
