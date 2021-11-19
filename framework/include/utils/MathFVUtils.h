//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "MooseError.h"
#include "FaceInfo.h"
#include "Limiter.h"
#include "MathUtils.h"
#include "libmesh/compare_types.h"
#include "libmesh/elem.h"
#include <tuple>

template <typename>
class MooseVariableFV;

namespace Moose
{
namespace FV
{
inline std::
    tuple<const FaceInfo *, Moose::FV::LimiterType, bool, bool, std::pair<SubdomainID, SubdomainID>>
    makeCDFace(const FaceInfo & fi,
               const std::pair<SubdomainID, SubdomainID> & subs,
               const bool skewness_correction = false)
{
  return std::make_tuple(
      &fi, Moose::FV::LimiterType::CentralDifference, true, skewness_correction, subs);
}

inline std::
    tuple<const FaceInfo *, Moose::FV::LimiterType, bool, bool, std::pair<SubdomainID, SubdomainID>>
    makeCDFace(const FaceInfo & fi, const bool skewness_correction = false)
{
  return makeCDFace(
      fi,
      std::make_pair(fi.elem().subdomain_id(),
                     fi.neighborPtr() ? fi.neighbor().subdomain_id() : Moose::INVALID_BLOCK_ID),
      skewness_correction);
}

/// This codifies a set of available ways to interpolate with elem+neighbor
/// solution information to calculate values (e.g. solution, material
/// properties, etc.) at the face (centroid).  These methods are used in the
/// class's interpolate functions.  Some interpolation methods are only meant
/// to be used with advective terms (e.g. upwind), others are more
/// generically applicable.
enum class InterpMethod
{
  /// (elem+neighbor)/2
  Average,
  /// (gc*elem+(1-gc)*neighbor)+gradient*(rf-rf')
  SkewCorrectedAverage,
  /// Extended stencil using the vertex values
  VertexBased,
  /// weighted
  Upwind,
  // Rhie-Chow
  RhieChow,
  VanLeer,
  MinMod,
  SOU,
  QUICK
};

/**
 * Produce the interpolation coefficients in the equation:
 *
 * \phi_f = c_1 * \phi_{F1} + c_2 * \phi_{F2}
 *
 * A couple of examples: if we are doing an average interpolation with
 * an orthogonal regular grid, then the pair will be (0.5, 0.5). If we are doing an
 * upwind interpolation with the velocity facing outward from the F1 element,
 * then the pair will be (1.0, 0.0).
 *
 * @param m The interpolation method
 * @param fi The face information
 * @param one_is_elem Whether fi.elem() == F1
 * @param advector The advecting velocity. Not relevant for an Average interpolation
 * @return a pair where the first Real is c_1 and the second Real is c_2
 */
template <typename Vector = RealVectorValue>
std::pair<Real, Real>
interpCoeffs(const InterpMethod m,
             const FaceInfo & fi,
             const bool one_is_elem,
             const Vector advector = Vector())
{
  switch (m)
  {
    case InterpMethod::Average:
    {
      if (one_is_elem)
        return std::make_pair(fi.gC(), 1. - fi.gC());
      else
        return std::make_pair(1. - fi.gC(), fi.gC());
    }

    case InterpMethod::SkewCorrectedAverage:
    {
      if (one_is_elem)
        return std::make_pair(fi.gCSkewed(), 1. - fi.gCSkewed());
      else
        return std::make_pair(1. - fi.gCSkewed(), fi.gCSkewed());
    }

    case InterpMethod::Upwind:
    {
      if ((advector * fi.normal() > 0) == one_is_elem)
        return std::make_pair(1., 0.);
      else
        return std::make_pair(0., 1.);
    }

    default:
      mooseError("Unrecognized interpolation method");
  }
}

/**
 * A simple linear interpolation of values between cell centers to a cell face. The \p one_is_elem
 * parameter indicates whether value1 corresponds to the FaceInfo elem value; else it corresponds to
 * the FaceInfo neighbor value
 */
template <typename T, typename T2>
typename libMesh::CompareTypes<T, T2>::supertype
linearInterpolation(const T & value1,
                    const T2 & value2,
                    const FaceInfo & fi,
                    const bool one_is_elem,
                    const InterpMethod interp_method = InterpMethod::Average)
{
  mooseAssert(interp_method == InterpMethod::Average ||
                  interp_method == InterpMethod::SkewCorrectedAverage,
              "The selected interpolation function only works with average or skewness-corrected "
              "average options!");
  const auto coeffs = interpCoeffs(interp_method, fi, one_is_elem);
  return coeffs.first * value1 + coeffs.second * value2;
}

/**
 * Linear interpolation with skewness correction using the face gradient.
 * See more info in Moukalled Chapter 9. The correction involves a first order
 * Taylor expansion around the intersection of the cell face and the line
 * connecting the two cell centers.
 */
template <typename T, typename T2, typename T3>
typename libMesh::CompareTypes<T, T2>::supertype
skewCorrectedlinearInterpolation(const T & value1,
                                 const T2 & value2,
                                 const T3 & face_gradient,
                                 const FaceInfo & fi,
                                 const bool one_is_elem)
{
  const auto coeffs = interpCoeffs(InterpMethod::SkewCorrectedAverage, fi, one_is_elem);

  auto value = (coeffs.first * value1 + coeffs.second * value2) +
               face_gradient * (fi.faceCentroid() - fi.rIntersection());
  return value;
}

/// Provides interpolation of face values for non-advection-specific purposes (although it can/will
/// still be used by advective kernels sometimes). The interpolated value is stored in result.
/// This should be called when a face value needs to be computed from two neighboring
/// cells/elements. value1 and value2 represent the cell property/values from which to compute the
/// face value. The \p one_is_elem parameter indicates whether value1 corresponds to the FaceInfo
/// elem value; else it corresponds to the FaceInfo neighbor value
template <typename T, typename T2, typename T3>
void
interpolate(InterpMethod m,
            T & result,
            const T2 & value1,
            const T3 & value2,
            const FaceInfo & fi,
            const bool one_is_elem)
{
  switch (m)
  {
    case InterpMethod::Average:
      result = linearInterpolation(value1, value2, fi, one_is_elem);
      break;
    default:
      mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
  }
}

/**
 * Computes the product of the advected and the advector based on the given interpolation method
 */
template <typename T1,
          typename T2,
          typename T3,
          template <typename>
          class Vector1,
          template <typename>
          class Vector2>
void
interpolate(InterpMethod m,
            Vector1<T1> & result,
            const T2 & fi_elem_advected,
            const T2 & fi_neighbor_advected,
            const Vector2<T3> & fi_elem_advector,
            const Vector2<T3> & fi_neighbor_advector,
            const FaceInfo & fi)
{
  switch (m)
  {
    case InterpMethod::Average:
      result = linearInterpolation(fi_elem_advected * fi_elem_advector,
                                   fi_neighbor_advected * fi_neighbor_advector,
                                   fi,
                                   true);
      break;
    case InterpMethod::Upwind:
    {
      const auto face_advector = linearInterpolation(MetaPhysicL::raw_value(fi_elem_advector),
                                                     MetaPhysicL::raw_value(fi_neighbor_advector),
                                                     fi,
                                                     true);
      Real elem_coeff, neighbor_coeff;
      if (face_advector * fi.normal() > 0)
        elem_coeff = 1, neighbor_coeff = 0;
      else
        elem_coeff = 0, neighbor_coeff = 1;

      result = elem_coeff * fi_elem_advected * fi_elem_advector +
               neighbor_coeff * fi_neighbor_advected * fi_neighbor_advector;
      break;
    }
    default:
      mooseError("unsupported interpolation method for FVFaceInterface::interpolate");
  }
}

/// Provides interpolation of face values for advective flux kernels.  This should be called by
/// advective kernels when a face value is needed from two neighboring cells/elements.  The
/// interpolated value is stored in result. value1 and value2 represent the two neighboring advected
/// cell property/values. advector represents the vector quantity at the face that is doing the
/// advecting (e.g. the flow velocity at the face); this value often will have been computed using a
/// call to the non-advective interpolate function. The \p one_is_elem parameter indicates whether
/// value1 corresponds to the FaceInfo elem value; else it corresponds to the FaceInfo neighbor
/// value
template <typename T, typename T2, typename T3, typename Vector>
void
interpolate(InterpMethod m,
            T & result,
            const T2 & value1,
            const T3 & value2,
            const Vector & advector,
            const FaceInfo & fi,
            const bool one_is_elem)
{
  const auto coeffs = interpCoeffs(m, fi, one_is_elem, advector);
  result = coeffs.first * value1 + coeffs.second * value2;
}

/// Calculates and returns "grad_u dot normal" on the face to be used for
/// diffusive terms.  If using any cross-diffusion corrections, etc. all
/// those calculations should be handled appropriately by this function.
template <typename T, typename T2>
ADReal gradUDotNormal(const T & elem_value,
                      const T2 & neighbor_value,
                      const FaceInfo & face_info,
                      const MooseVariableFV<Real> & fv_var,
                      bool correct_skewness = false);

/**
 * From Moukalled 12.30
 *
 * r_f = (phiC - phiU) / (phiD - phiC)
 *
 * However, this formula is only clear on grids where upwind locations can be readily determined,
 * which is not the case for unstructured meshes. So we leverage a virtual upwind location and
 * Moukalled 12.65
 *
 * phiD - phiU = 2 * grad(phi)_C * dCD ->
 * phiU = phiD - 2 * grad(phi)_C * dCD
 *
 * Combining the two equations and performing some algebraic manipulation yields this equation for
 * r_f:
 *
 * r_f = 2 * grad(phi)_C * dCD / (phiD - phiC) - 1
 *
 * This equation is clearly asymmetric considering the face between C and D because of the
 * subscript on grad(phi). Hence this method can be thought of as constructing an r associated with
 * the C side of the face
 */
template <typename Scalar, typename Vector>
Scalar
rF(const Scalar & phiC, const Scalar & phiD, const Vector & gradC, const RealVectorValue & dCD)
{
  if ((phiD - phiC) == 0)
    return 1e6 * MathUtils::sign(gradC * dCD) + 0 * (gradC * dCD + phiD + phiC);

  return 2. * gradC * dCD / (phiD - phiC) - 1.;
}

/**
 * Produce the interpolation coefficients in the equation:
 *
 * \phi_f = c_upwind * \phi_{upwind} + c_downwind * \phi_{downwind}
 *
 * A couple of examples: if we are doing an average interpolation with
 * an orthogonal regular grid, then the pair will be (0.5, 0.5). If we are doing an
 * upwind interpolation then the pair will be (1.0, 0.0).
 *
 * @return a pair where the first Real is c_upwind and the second Real is c_downwind
 */
template <typename T>
std::pair<T, T>
interpCoeffs(const Limiter<T> & limiter,
             const T & phi_upwind,
             const T & phi_downwind,
             const VectorValue<T> * const grad_phi_upwind,
             const FaceInfo & fi,
             const bool fi_elem_is_upwind)
{
  // Using beta, w_f, g nomenclature from Greenshields
  const auto beta = limiter(phi_upwind,
                            phi_downwind,
                            grad_phi_upwind,
                            fi_elem_is_upwind ? fi.dCF() : RealVectorValue(-fi.dCF()));

  const auto w_f = fi_elem_is_upwind ? fi.gC() : (1. - fi.gC());

  const auto g = beta * (1. - w_f);

  return std::make_pair(1. - g, g);
}

template <typename T>
struct GradientType
{
};

/**
 * Interpolates with a limiter
 */
template <typename Scalar,
          typename Vector,
          typename Enable = typename std::enable_if<ScalarTraits<Scalar>::value>::type>
Scalar
interpolate(const Limiter<Scalar> & limiter,
            const Scalar & phi_upwind,
            const Scalar & phi_downwind,
            const Vector * const grad_phi_upwind,
            const FaceInfo & fi,
            const bool fi_elem_is_upwind)
{
  auto pr = interpCoeffs(limiter, phi_upwind, phi_downwind, grad_phi_upwind, fi, fi_elem_is_upwind);
  return pr.first * phi_upwind + pr.second * phi_downwind;
}

/**
 * Vector overload
 */
template <typename Limiter, typename T, typename Tensor>
VectorValue<T>
interpolate(const Limiter & limiter,
            const TypeVector<T> & phi_upwind,
            const TypeVector<T> & phi_downwind,
            const Tensor * const /*grad_phi_upwind*/,
            const FaceInfo & fi,
            const bool fi_elem_is_upwind)
{
  mooseAssert(limiter.constant(),
              "Non-constant limiters are not currently supported in the vector overload of the "
              "limited interpolate method");
  static const VectorValue<T> example_gradient(0);

  VectorValue<T> ret;
  for (const auto i : make_range(unsigned(LIBMESH_DIM)))
    ret(i) = interpolate(
        limiter, phi_upwind(i), phi_downwind(i), &example_gradient, fi, fi_elem_is_upwind);

  return ret;
}
}
}
