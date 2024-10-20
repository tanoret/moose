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
#include "MooseFunctor.h"
#include "libmesh/compare_types.h"
#include "libmesh/elem.h"
#include <tuple>

template <typename>
class MooseVariableFV;
class FaceArgInterface;

namespace Moose
{
namespace FV
{
/// This codifies a set of available ways to interpolate with elem+neighbor
/// solution information to calculate values (e.g. solution, material
/// properties, etc.) at the face (centroid).  These methods are used in the
/// class's interpolate functions.  Some interpolation methods are only meant
/// to be used with advective terms (e.g. upwind), others are more
/// generically applicable.
enum class InterpMethod
{
  /// gc*elem+(1-gc)*neighbor
  Average,
  /// 1/(gc/elem+(1-gc)/neighbor)
  HarmonicAverage,
  /// (gc*elem+(1-gc)*neighbor)+gradient*(rf-rf')
  SkewCorrectedAverage,
  /// weighted
  Upwind,
  // Rhie-Chow
  RhieChow,
  VanLeer,
  MinMod,
  SOU,
  QUICK,
  Venkatakrishnan
};

enum class LinearFVComputationMode
{
  RHS,
  Matrix,
  FullSystem
};

/**
 * Returns an enum with all the currently supported interpolation methods and the current default
 * for FV: first-order upwind
 * @return MooseEnum with all the face interpolation methods supported
 */
MooseEnum interpolationMethods();

/**
 * @return An input parameters object that contains the \p advected_interp_method parameter, e.g.
 * the interpolation method to use for an advected quantity
 */
InputParameters advectedInterpolationParameter();

/*
 * Converts from the interpolation method to the interpolation enum.
 * This routine is here in lieu of using a MooseEnum for InterpMethod
 * @param interp_method the name of the interpolation method
 * @return the interpolation method
 */
InterpMethod selectInterpolationMethod(const std::string & interp_method);

/**
 * Sets one interpolation method
 * @param obj The \p MooseObject with input parameters to query
 * @param interp_method The interpolation method we will set
 * @param param_name The name of the parameter setting this interpolation method
 * @return Whether the interpolation method has indicated that we will need more than the
 * default level of ghosting
 */
bool setInterpolationMethod(const MooseObject & obj,
                            Moose::FV::InterpMethod & interp_method,
                            const std::string & param_name);

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
 * The template is needed in case the face flux carries derivatives in an AD setting.
 *
 * @param m The interpolation method
 * @param fi The face information
 * @param one_is_elem Whether fi.elem() == F1
 * @param face_flux The advecting face flux. Not relevant for an Average interpolation
 * @return a pair where the first Real is c_1 and the second Real is c_2
 */
template <typename T = Real>
std::pair<Real, Real>
interpCoeffs(const InterpMethod m,
             const FaceInfo & fi,
             const bool one_is_elem,
             const T & face_flux = 0.0)
{
  switch (m)
  {
    case InterpMethod::Average:
    case InterpMethod::SkewCorrectedAverage:
    {
      if (one_is_elem)
        return std::make_pair(fi.gC(), 1. - fi.gC());
      else
        return std::make_pair(1. - fi.gC(), fi.gC());
    }

    case InterpMethod::Upwind:
    {
      if ((face_flux > 0) == one_is_elem)
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
 * Computes the harmonic mean (1/(gc/value1+(1-gc)/value2)) of Reals, RealVectorValues and
 * RealTensorValues while accounting for the possibility that one or both of them are AD.
 * For tensors, we use a component-wise mean instead of the matrix-inverse based option.
 * @param value1 Reference to value1 in the (1/(gc/value1+(1-gc)/value2)) expression
 * @param value2 Reference to value2 in the (1/(gc/value1+(1-gc)/value2)) expression
 * @param fi Reference to the FaceInfo of the face on which the interpolation is requested
 * @param one_is_elem Boolean indicating if the interpolation weight on FaceInfo belongs to the
 * element corresponding to value1
 * @return The interpolated value
 */
template <typename T1, typename T2>
typename libMesh::CompareTypes<T1, T2>::supertype
harmonicInterpolation(const T1 & value1,
                      const T2 & value2,
                      const FaceInfo & fi,
                      const bool one_is_elem)
{
  // We check if the base values of the given template types match, if not we throw a compile-time
  // error
  static_assert(std::is_same<typename MetaPhysicL::RawType<T1>::value_type,
                             typename MetaPhysicL::RawType<T2>::value_type>::value,
                "The input values for harmonic interpolation need to have the same raw-value!");

  // Fetch the interpolation coefficients, we use the exact same coefficients as for a simple
  // weighted average
  const auto coeffs = interpCoeffs(InterpMethod::Average, fi, one_is_elem);

  // We check if the types are fit to compute the harmonic mean of. This is done compile-time
  // using constexpr. We start with Real/ADReal which is straightforward if the input values are
  // positive.
  if constexpr (libMesh::TensorTools::TensorTraits<T1>::rank == 0)
  {
    // The harmonic mean of mixed positive and negative numbers (and 0 as well) is not well-defined
    // so we assert that the input values shall be positive.
#ifndef NDEBUG
    if (value1 * value2 <= 0)
      mooseWarning("Input values (" + Moose::stringify(MetaPhysicL::raw_value(value1)) + " & " +
                   Moose::stringify(MetaPhysicL::raw_value(value2)) +
                   ") must be of the same sign for harmonic interpolation");
#endif
    return 1.0 / (coeffs.first / value1 + coeffs.second / value2);
  }
  // For vectors (ADRealVectorValue, VectorValue), we take the component-wise harmonic mean
  else if constexpr (libMesh::TensorTools::TensorTraits<T1>::rank == 1)
  {
    typename libMesh::CompareTypes<T1, T2>::supertype result;
    for (const auto i : make_range(Moose::dim))
    {
#ifndef NDEBUG
      if (value1(i) * value2(i) <= 0)
        mooseWarning("Component " + std::to_string(i) + " of input values (" +
                     Moose::stringify(MetaPhysicL::raw_value(value1(i))) + " & " +
                     Moose::stringify(MetaPhysicL::raw_value(value2(i))) +
                     ") must be of the same sign for harmonic interpolation");
#endif
      result(i) = 1.0 / (coeffs.first / value1(i) + coeffs.second / value2(i));
    }
    return result;
  }
  // For tensors (ADRealTensorValue, TensorValue), similarly to the vectors,
  // we take the component-wise harmonic mean instead of the matrix-inverse approach
  else if constexpr (libMesh::TensorTools::TensorTraits<T1>::rank == 2)
  {
    typename libMesh::CompareTypes<T1, T2>::supertype result;
    for (const auto i : make_range(Moose::dim))
      for (const auto j : make_range(Moose::dim))
      {
#ifndef NDEBUG
        if (value1(i, j) * value2(i, j) <= 0)
          mooseWarning("Component (" + std::to_string(i) + "," + std::to_string(j) +
                       ") of input values (" +
                       Moose::stringify(MetaPhysicL::raw_value(value1(i, j))) + " & " +
                       Moose::stringify(MetaPhysicL::raw_value(value2(i, j))) +
                       ") must be of the same sign for harmonic interpolation");
#endif
        result(i, j) = 1.0 / (coeffs.first / value1(i, j) + coeffs.second / value2(i, j));
      }
    return result;
  }
  // We ran out of options, harmonic mean is not supported for other types at the moment
  else
    // This line is supposed to throw an error when the user tries to compile this function with
    // types that are not supported. This is the reason we needed the always_false function. Hope as
    // C++ gets nicer, we can do this in a nicer way.
    static_assert(Moose::always_false<T1>,
                  "Harmonic interpolation is not implemented for the used type!");
}

/**
 * Linear interpolation with skewness correction using the face gradient.
 * See more info in Moukalled Chapter 9. The correction involves a first order
 * Taylor expansion around the intersection of the cell face and the line
 * connecting the two cell centers.
 */
template <typename T, typename T2, typename T3>
typename libMesh::CompareTypes<T, T2>::supertype
skewCorrectedLinearInterpolation(const T & value1,
                                 const T2 & value2,
                                 const T3 & face_gradient,
                                 const FaceInfo & fi,
                                 const bool one_is_elem)
{
  const auto coeffs = interpCoeffs(InterpMethod::SkewCorrectedAverage, fi, one_is_elem);

  auto value = (coeffs.first * value1 + coeffs.second * value2) +
               face_gradient * fi.skewnessCorrectionVector();
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
    case InterpMethod::SkewCorrectedAverage:
      result = linearInterpolation(value1, value2, fi, one_is_elem);
      break;
    case InterpMethod::HarmonicAverage:
      result = harmonicInterpolation(value1, value2, fi, one_is_elem);
      break;
    default:
      mooseError("unsupported interpolation method for interpolate() function");
  }
}

/**
 * perform a possibly skew-corrected linear interpolation by evaluating the supplied functor with
 * the provided functor face argument
 */
template <typename T, FunctorEvaluationKind FEK = FunctorEvaluationKind::Value>
T
linearInterpolation(const FunctorBase<T> & functor, const FaceArg & face, const StateArg & time)
{
  static_assert((FEK == FunctorEvaluationKind::Value) || (FEK == FunctorEvaluationKind::Dot),
                "Only Value and Dot evaluations currently supported");
  mooseAssert(face.limiter_type == LimiterType::CentralDifference,
              "this interpolation method is meant for linear interpolations");

  mooseAssert(face.fi,
              "We must have a non-null face_info in order to prepare our ElemFromFace tuples");

  constexpr FunctorEvaluationKind GFEK = FunctorGradientEvaluationKind<FEK>::value;

  const auto elem_arg = face.makeElem();
  const auto neighbor_arg = face.makeNeighbor();

  const auto elem_eval = functor.template genericEvaluate<FEK>(elem_arg, time);
  const auto neighbor_eval = functor.template genericEvaluate<FEK>(neighbor_arg, time);

  if (face.correct_skewness)
  {
    // This condition ensures that the recursive algorithm (face_center->
    // face_gradient -> cell_gradient -> face_center -> ...) terminates after
    // one loop. It is hardcoded to one loop at this point since it yields
    // 2nd order accuracy on skewed meshes with the minimum additional effort.
    FaceArg new_face(face);
    new_face.correct_skewness = false;
    const auto surface_gradient = functor.template genericEvaluate<GFEK>(new_face, time);

    return skewCorrectedLinearInterpolation(
        elem_eval, neighbor_eval, surface_gradient, *face.fi, true);
  }
  else
    return linearInterpolation(elem_eval, neighbor_eval, *face.fi, true);
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
  const auto coeffs = interpCoeffs(m, fi, one_is_elem, advector * fi.normal());
  result = coeffs.first * value1 + coeffs.second * value2;
}

/// Calculates and returns "grad_u dot normal" on the face to be used for
/// diffusive terms.  If using any cross-diffusion corrections, etc. all
/// those calculations should be handled appropriately by this function.
ADReal gradUDotNormal(const FaceInfo & face_info,
                      const MooseVariableFV<Real> & fv_var,
                      const Moose::StateArg & time,
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
  static const auto zero_vec = RealVectorValue(0);
  if ((phiD - phiC) == 0)
    // Handle zero denominator case. Note that MathUtils::sign returns 1 for sign(0) so we can omit
    // that operation here (e.g. sign(phiD - phiC) = sign(0) = 1). The second term preserves the
    // same sparsity pattern as the else branch; we want to add this so that we don't risk PETSc
    // shrinking the matrix now and then potentially reallocating nonzeros later (which is very
    // slow)
    return 1e6 * MathUtils::sign(gradC * dCD) + zero_vec * gradC;

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
  const auto beta = limiter(
      phi_upwind, phi_downwind, grad_phi_upwind, fi_elem_is_upwind ? fi.dCN() : Point(-fi.dCN()));

  const auto w_f = fi_elem_is_upwind ? fi.gC() : (1. - fi.gC());

  const auto g = beta * (1. - w_f);

  return std::make_pair(1. - g, g);
}

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
            const Tensor * const grad_phi_upwind,
            const FaceInfo & fi,
            const bool fi_elem_is_upwind)
{
  mooseAssert(limiter.constant() || grad_phi_upwind,
              "Non-null gradient only supported for constant limiters.");

  const VectorValue<T> * const gradient_example = nullptr;
  VectorValue<T> ret;
  for (const auto i : make_range(unsigned(LIBMESH_DIM)))
  {
    if (grad_phi_upwind)
    {
      const VectorValue<T> gradient = grad_phi_upwind->row(i);
      ret(i) =
          interpolate(limiter, phi_upwind(i), phi_downwind(i), &gradient, fi, fi_elem_is_upwind);
    }
    else
      ret(i) = interpolate(
          limiter, phi_upwind(i), phi_downwind(i), gradient_example, fi, fi_elem_is_upwind);
  }

  return ret;
}

/**
 * @brief Computes the full limited interpolation of a variable using a specified limiter.
 *
 * This function performs a full limited interpolation of a variable, taking into account
 * the values and gradients at both upwind and downwind locations, as well as geometric
 * information and limits. It applies the specified limiter to ensure the interpolation
 * adheres to the constraints imposed by the limiter.
 *
 * @tparam T The data type of the scalar values and the return type.
 * @param limiter The limiter object used to compute the flux limiting ratio.
 * @param phi_upwind The scalar value at the upwind location.
 * @param phi_downwind The scalar value at the downwind location.
 * @param grad_phi_upwind Pointer to the gradient vector at the upwind location.
 * @param grad_phi_downwind Pointer to the gradient vector at the downwind location.
 * @param fi FaceInfo object containing geometric details such as face centroid and cell centroids.
 * @param fi_elem_is_upwind Boolean indicating if the face info element is upwind.
 * @param max_value The maximum allowable value.
 * @param min_value The minimum allowable value.
 * @return The computed limited interpolation value.
 *
 * This function performs the following steps:
 * 1. Computes the direction vector \( dCD \) based on whether the current element is upwind.
 * 2. Calls the `limiter` functor with the provided scalar values, gradient vectors, direction
 * vector, and geometric information to compute the flux limiting ratio \( \beta \).
 * 3. Determines the face centroid and the appropriate cell centroid based on whether the current
 * element is upwind.
 * 4. Computes the delta value at the face by taking the dot product of the upwind gradient vector
 * with the difference between the face and cell centroids.
 * 5. Returns the interpolated value by adding the product of the limiter \( \beta \) and the delta
 * value to the upwind scalar value.
 */
template <typename T>
T
fullLimitedInterpolation(const Limiter<T> & limiter,
                         const T & phi_upwind,
                         const T & phi_downwind,
                         const VectorValue<T> * const grad_phi_upwind,
                         const VectorValue<T> * const grad_phi_downwind,
                         const FaceInfo & fi,
                         const bool fi_elem_is_upwind,
                         const T & max_value,
                         const T & min_value)
{
  // Compute the direction vector based on whether the current element is upwind
  const auto dCD = fi_elem_is_upwind ? fi.dCN() : Point(-fi.dCN());

  // Compute the flux limiting ratio using the specified limiter
  const auto beta = limiter(phi_upwind,
                            phi_downwind,
                            grad_phi_upwind,
                            grad_phi_downwind,
                            dCD,
                            max_value,
                            min_value,
                            &fi,
                            fi_elem_is_upwind);

  // Determine the face centroid and the appropriate cell centroid
  const auto & face_centroid = fi.faceCentroid();
  const auto & cell_centroid = fi_elem_is_upwind ? fi.elemCentroid() : fi.neighborCentroid();

  // Compute the delta value at the face
  const auto & delta_face = (*grad_phi_upwind) * (face_centroid - cell_centroid);

  // Return the interpolated value
  return phi_upwind + beta * delta_face;
}

/**
 * Interpolates with a limiter and a face argument
 * @return a pair of pairs. The first pair corresponds to the interpolation coefficients with the
 * first corresponding to the face information element and the second corresponding to the face
 * information neighbor. This pair should sum to unity. The second pair corresponds to the face
 * information functor element value and neighbor
 */
template <typename T,
          FunctorEvaluationKind FEK = FunctorEvaluationKind::Value,
          typename Enable = typename std::enable_if<ScalarTraits<T>::value>::type>
std::pair<std::pair<T, T>, std::pair<T, T>>
interpCoeffsAndAdvected(const FunctorBase<T> & functor, const FaceArg & face, const StateArg & time)
{
  static_assert((FEK == FunctorEvaluationKind::Value) || (FEK == FunctorEvaluationKind::Dot),
                "Only Value and Dot evaluations currently supported");

  constexpr FunctorEvaluationKind GFEK = FunctorGradientEvaluationKind<FEK>::value;
  typedef typename FunctorBase<T>::GradientType GradientType;
  static const GradientType zero(0);

  mooseAssert(face.fi, "this must be non-null");
  const auto limiter = Limiter<typename LimiterValueType<T>::value_type>::build(face.limiter_type);

  const auto upwind_arg = face.elem_is_upwind ? face.makeElem() : face.makeNeighbor();
  const auto downwind_arg = face.elem_is_upwind ? face.makeNeighbor() : face.makeElem();
  auto phi_upwind = functor.template genericEvaluate<FEK>(upwind_arg, time);
  auto phi_downwind = functor.template genericEvaluate<FEK>(downwind_arg, time);

  std::pair<T, T> interp_coeffs;
  if (face.limiter_type == LimiterType::Upwind ||
      face.limiter_type == LimiterType::CentralDifference)
    interp_coeffs =
        interpCoeffs(*limiter, phi_upwind, phi_downwind, &zero, *face.fi, face.elem_is_upwind);
  else
  {
    const auto grad_phi_upwind = functor.template genericEvaluate<GFEK>(upwind_arg, time);
    interp_coeffs = interpCoeffs(
        *limiter, phi_upwind, phi_downwind, &grad_phi_upwind, *face.fi, face.elem_is_upwind);
  }

  if (face.elem_is_upwind)
    return std::make_pair(std::move(interp_coeffs),
                          std::make_pair(std::move(phi_upwind), std::move(phi_downwind)));
  else
    return std::make_pair(
        std::make_pair(std::move(interp_coeffs.second), std::move(interp_coeffs.first)),
        std::make_pair(std::move(phi_downwind), std::move(phi_upwind)));
}

// Utility function to compute the min-max values in a two-cell stencil
template <typename T,
          FunctorEvaluationKind FEK = FunctorEvaluationKind::Value,
          typename Enable = typename std::enable_if<ScalarTraits<T>::value>::type>
std::pair<T, T>
computeMinMaxValue(const FunctorBase<T> & functor, const FaceArg & face, const StateArg & time)
{
  T max_value(0), min_value(1e10);
  for (auto neighbor : (*face.fi).elem().neighbor_ptr_range())
  {
    // If not a valid neighbor
    if (neighbor == nullptr)
      continue;
    else
    {
      const ElemArg local_cell_arg = {neighbor, false};
      const auto value = functor.template genericEvaluate<FEK>(local_cell_arg, time);
      max_value = std::max(max_value, value);
      min_value = std::min(min_value, value);
    }
  }
  for (auto neighbor : (*face.fi).neighbor().neighbor_ptr_range())
  {
    // If not a valid neighbor
    if (neighbor == nullptr)
      continue;
    else
    {
      const ElemArg local_cell_arg = {neighbor, false};
      const auto value = functor.template genericEvaluate<FEK>(local_cell_arg, time);
      max_value = std::max(max_value, value);
      min_value = std::min(min_value, value);
    }
  }
  return std::make_pair(max_value, min_value);
}

template <typename T,
          FunctorEvaluationKind FEK = FunctorEvaluationKind::Value,
          typename Enable = typename std::enable_if<ScalarTraits<T>::value>::type>
T
interpolate(const FunctorBase<T> & functor, const FaceArg & face, const StateArg & time)
{
  static_assert((FEK == FunctorEvaluationKind::Value) || (FEK == FunctorEvaluationKind::Dot),
                "Only Value and Dot evaluations currently supported");

  // Special handling for central differencing as it is the only interpolation method which
  // currently supports skew correction
  if (face.limiter_type == LimiterType::CentralDifference)
    return linearInterpolation<T, FEK>(functor, face, time);

  if (face.limiter_type == LimiterType::Upwind ||
      face.limiter_type == LimiterType::CentralDifference)
  {
    const auto [interp_coeffs, advected] = interpCoeffsAndAdvected<T, FEK>(functor, face, time);
    return interp_coeffs.first * advected.first + interp_coeffs.second * advected.second;
  }
  else
  {

    constexpr FunctorEvaluationKind GFEK = FunctorGradientEvaluationKind<FEK>::value;
    typedef typename FunctorBase<T>::GradientType GradientType;
    static const GradientType zero(0);
    const auto & limiter =
        Limiter<typename LimiterValueType<T>::value_type>::build(face.limiter_type);

    const auto & upwind_arg = face.elem_is_upwind ? face.makeElem() : face.makeNeighbor();
    const auto & downwind_arg = face.elem_is_upwind ? face.makeNeighbor() : face.makeElem();
    const auto & phi_upwind = functor.template genericEvaluate<FEK>(upwind_arg, time);
    const auto & phi_downwind = functor.template genericEvaluate<FEK>(downwind_arg, time);

    const auto & time_arg = Moose::StateArg(1, Moose::SolutionIterationType::Time);

    // Min-Max value
    T max_value, min_value;
    std::tie(max_value, min_value) = computeMinMaxValue(functor, face, time_arg);

    const auto & grad_phi_upwind = functor.template genericEvaluate<GFEK>(upwind_arg, time_arg);
    const auto & grad_phi_downwind = functor.template genericEvaluate<GFEK>(face, time_arg);

    return fullLimitedInterpolation(*limiter,
                                    phi_upwind,
                                    phi_downwind,
                                    &grad_phi_upwind,
                                    &grad_phi_downwind,
                                    *face.fi,
                                    face.elem_is_upwind,
                                    max_value,
                                    min_value);
  }
}

// Utility function to compute the min-max values in a two-cell stencil
template <typename T>
std::pair<T, T>
computeMinMaxValue(const FunctorBase<VectorValue<T>> & functor,
                   const FaceArg & face,
                   const StateArg & time,
                   const unsigned int & component)
{
  T max_value(0), min_value(1e10);
  for (auto neighbor : (*face.fi).elem().neighbor_ptr_range())
  {
    // If not a valid neighbor
    if (neighbor == nullptr)
      continue;
    else
    {
      const ElemArg local_cell_arg = {neighbor, false};
      const auto value = functor(local_cell_arg, time)(component);
      max_value = std::max(max_value, value);
      min_value = std::min(min_value, value);
    }
  }
  for (auto neighbor : (*face.fi).neighbor().neighbor_ptr_range())
  {
    // If not a valid neighbor
    if (neighbor == nullptr)
      continue;
    else
    {
      const ElemArg local_cell_arg = {neighbor, false};
      const auto value = functor(local_cell_arg, time)(component);
      max_value = std::max(max_value, value);
      min_value = std::min(min_value, value);
    }
  }
  return std::make_pair(max_value, min_value);
}

template <typename T>
VectorValue<T>
interpolate(const FunctorBase<VectorValue<T>> & functor,
            const FaceArg & face,
            const StateArg & time)
{
  static const VectorValue<T> grad_zero(0);

  mooseAssert(face.fi, "this must be non-null");
  const auto limiter = Limiter<typename LimiterValueType<T>::value_type>::build(face.limiter_type);

  const auto upwind_arg = face.elem_is_upwind ? face.makeElem() : face.makeNeighbor();
  const auto downwind_arg = face.elem_is_upwind ? face.makeNeighbor() : face.makeElem();
  auto phi_upwind = functor(upwind_arg, time);
  auto phi_downwind = functor(downwind_arg, time);

  VectorValue<T> ret;
  T coeff_upwind, coeff_downwind;

  if (face.limiter_type == LimiterType::Upwind ||
      face.limiter_type == LimiterType::CentralDifference)
  {
    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      const auto &component_upwind = phi_upwind(i), component_downwind = phi_downwind(i);
      std::tie(coeff_upwind, coeff_downwind) = interpCoeffs(*limiter,
                                                            component_upwind,
                                                            component_downwind,
                                                            &grad_zero,
                                                            *face.fi,
                                                            face.elem_is_upwind);
      ret(i) = coeff_upwind * component_upwind + coeff_downwind * component_downwind;
    }
  }
  else
  {
    const auto & time_arg = Moose::StateArg(1, Moose::SolutionIterationType::Time);
    const auto & grad_phi_upwind = functor.gradient(upwind_arg, time_arg);
    const auto & grad_phi_downwind = functor.gradient(downwind_arg, time_arg);

    for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
    {
      const auto &component_upwind = phi_upwind(i), component_downwind = phi_downwind(i);
      const VectorValue<T> &grad_upwind = grad_phi_upwind.row(i),
                           grad_downwind = grad_phi_downwind.row(i);

      // Min-Max value
      T max_value, min_value;
      std::tie(max_value, min_value) = computeMinMaxValue(functor, face, time_arg, i);

      ret(i) = fullLimitedInterpolation(*limiter,
                                        component_upwind,
                                        component_downwind,
                                        &grad_upwind,
                                        &grad_downwind,
                                        *face.fi,
                                        face.elem_is_upwind,
                                        max_value,
                                        min_value);
    }
  }

  return ret;
}

// Utility function to compute the min-max values in a two-cell stencil
template <typename T>
std::pair<T, T>
containerComputeMinMaxValue(const FunctorBase<T> & functor,
                            const FaceArg & face,
                            const StateArg & time,
                            const unsigned int & component)
{
  T max_value(0), min_value(1e10);
  for (auto neighbor : (*face.fi).elem().neighbor_ptr_range())
  {
    // If not a valid neighbor
    if (neighbor == nullptr)
      continue;
    else
    {
      const ElemArg local_cell_arg = {neighbor, false};
      const auto value = functor(local_cell_arg, time)[component];
      max_value = std::max(max_value, value);
      min_value = std::min(min_value, value);
    }
  }
  for (auto neighbor : (*face.fi).neighbor().neighbor_ptr_range())
  {
    // If not a valid neighbor
    if (neighbor == nullptr)
      continue;
    else
    {
      const ElemArg local_cell_arg = {neighbor, false};
      const auto value = functor(local_cell_arg, time)[component];
      max_value = std::max(max_value, value);
      min_value = std::min(min_value, value);
    }
  }
  return std::make_pair(max_value, min_value);
}

template <typename T>
T
containerInterpolate(const FunctorBase<T> & functor, const FaceArg & face, const StateArg & time)
{
  typedef typename FunctorBase<T>::GradientType ContainerGradientType;
  typedef typename ContainerGradientType::value_type GradientType;
  const GradientType * const example_gradient = nullptr;

  mooseAssert(face.fi, "this must be non-null");
  const auto limiter = Limiter<typename T::value_type>::build(face.limiter_type);

  const auto upwind_arg = face.elem_is_upwind ? face.makeElem() : face.makeNeighbor();
  const auto downwind_arg = face.elem_is_upwind ? face.makeNeighbor() : face.makeElem();
  const auto phi_upwind = functor(upwind_arg, time);
  const auto phi_downwind = functor(downwind_arg, time);

  // initialize in order to get proper size
  T ret = phi_upwind;
  typename T::value_type coeff_upwind, coeff_downwind;

  if (face.limiter_type == LimiterType::Upwind ||
      face.limiter_type == LimiterType::CentralDifference)
  {
    for (const auto i : index_range(ret))
    {
      const auto &component_upwind = phi_upwind[i], component_downwind = phi_downwind[i];
      std::tie(coeff_upwind, coeff_downwind) = interpCoeffs(*limiter,
                                                            component_upwind,
                                                            component_downwind,
                                                            example_gradient,
                                                            *face.fi,
                                                            face.elem_is_upwind);
      ret[i] = coeff_upwind * component_upwind + coeff_downwind * component_downwind;
    }
  }
  else
  {
    const auto & time_arg = Moose::StateArg(1, Moose::SolutionIterationType::Time);
    const auto & grad_phi_upwind = functor.gradient(upwind_arg, time_arg);
    const auto & grad_phi_downwind = functor.gradient(downwind_arg, time_arg);

    for (const auto i : index_range(ret))
    {
      const auto &component_upwind = phi_upwind[i], component_downwind = phi_downwind[i];
      const VectorValue<T> &grad_upwind = grad_phi_upwind[i], grad_downwind = grad_phi_downwind[i];

      // Min-Max value
      T max_value, min_value;
      std::tie(max_value, min_value) = containerComputeMinMaxValue(functor, face, time_arg, i);

      ret[i] = fullLimitedInterpolation(*limiter,
                                        component_upwind,
                                        component_downwind,
                                        &grad_upwind,
                                        &grad_downwind,
                                        *face.fi,
                                        face.elem_is_upwind,
                                        max_value,
                                        min_value);
    }
  }

  return ret;
}

template <typename T>
std::vector<T>
interpolate(const FunctorBase<std::vector<T>> & functor,
            const FaceArg & face,
            const StateArg & time)
{
  return containerInterpolate(functor, face, time);
}

template <typename T, std::size_t N>
std::array<T, N>
interpolate(const FunctorBase<std::array<T, N>> & functor,
            const FaceArg & face,
            const StateArg & time)
{
  return containerInterpolate(functor, face, time);
}

/**
 * Return whether the supplied face is on a boundary of the \p object's execution
 */
template <typename SubdomainRestrictable>
bool
onBoundary(const SubdomainRestrictable & obj, const FaceInfo & fi)
{
  const bool internal = fi.neighborPtr() && obj.hasBlocks(fi.elemSubdomainID()) &&
                        obj.hasBlocks(fi.neighborSubdomainID());
  return !internal;
}

/**
 * Determine whether the passed-in face is on the boundary of an object that lives on the provided
 * subdomains. Note that if the subdomain set is empty we consider that to mean that the object has
 * no block restriction and lives on the entire mesh
 */
bool onBoundary(const std::set<SubdomainID> & subs, const FaceInfo & fi);
}
}
