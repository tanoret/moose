//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <tuple>

#include "MooseMesh.h"
#include "MooseTypes.h"
#include "MooseError.h"
#include "Limiter.h"

#include "libmesh/elem.h"
#include "libmesh/quadrature.h"
#include "libmesh/remote_elem.h"
#include "libmesh/tensor_tools.h"

#include <unordered_map>
#include <functional>

class FaceInfo;

namespace Moose
{
/**
 * A non-templated base class for functors that allow an owner object to hold
 * different class template instantiations of \p Functor in a single container
 */
class FunctorBase
{
public:
  FunctorBase() = default;
  virtual ~FunctorBase() = default;

  ///@{
  /**
   * Virtual methods meant to be used for handling functor evaluation cache clearance
   */
  virtual void timestepSetup() = 0;
  virtual void residualSetup() = 0;
  virtual void jacobianSetup() = 0;
  ///@}
};

/**
 * Base class template for functor objects. This class template defines various \p operator()
 * overloads that allow a user to evaluate the functor at arbitrary geometric locations. This
 * template is meant to enable highly flexible on-the-fly variable and material property evaluations
 */
template <typename T>
class FunctorImpl : public FunctorBase
{
public:
  /**
   * A typedef defining a "face" evaluation calling argument. This is composed of
   * - a face information object which defines our location in space
   * - a limiter which defines how the functor evaluated on either side of the face should be\n
   *   interpolated to the face
   * - a boolean which states whether the face information element is upwind of the face
   * - a boolean indicating whether to perform skew correction
   * - a pair of subdomain IDs. These do not always correspond to the face info element subdomain\n
   *   ID and face info neighbor subdomain ID. For instance if a flux kernel is operating at a\n
   *   subdomain boundary on which the kernel is defined on one side but not the other, the\n
   *   passed-in subdomain IDs will both correspond to the subdomain ID that the flux kernel is\n
   *   defined on
   */
  using FaceArg = std::tuple<const FaceInfo *,
                             Moose::FV::LimiterType,
                             bool,
                             bool,
                             std::pair<SubdomainID, SubdomainID>>;

  /**
   * A typedef defining a "single-sided face" evaluation calling argument. This is identical to the
   * \p FaceArg argument with the exception that a single SubdomainID is given as the last argument
   * to the tuple. This allows disambiguation, in \p FunctorMaterialProperty contexts, of which
   * property definition to use at a face on which there may be different property definitions on
   * either side of the face. Additionally this overload can be useful when on an external boundary
   * face when we want to avoid using ghost information
   */
  using SingleSidedFaceArg =
      std::tuple<const FaceInfo *, Moose::FV::LimiterType, bool, SubdomainID>;

  /**
   * People should think of this geometric argument as corresponding to the location in space of the
   * provided element centroid, \b not as corresponding to the location of the provided face
   * information. Summary of data in this argument:
   * - an element, whose centroid we should think of as the evaluation point. It is possible that\n
   *   the element will be a nullptr in which case, the evaluation point should be thought of as\n
   *   the location of a ghosted element centroid
   * - a face information object. When the provided element is null or for instance when the\n
   *   functoris a variable that does not exist on the provided element subdomain, this face\n
   *   information object will be used to help construct a ghost value evaluation
   * - a subdomain ID. This is useful when the functor is a material property and the user wants\n
   *   to indicate which material property definition should be used to evaluate the functor. For\n
   *   instance if we are using a flux kernel that is not defined on one side of the face, the\n
   *   subdomain ID will allow us to compute a ghost material property evaluation
   */
  using ElemFromFaceArg = std::tuple<const libMesh::Elem *, const FaceInfo *, SubdomainID>;

  /**
   * Argument for requesting functor evaluation at a quadrature point location in an element. Data
   * in the argument:
   * - The element containing the quadrature point
   * - The quadrature point index, e.g. if there are \p n quadrature points, we are requesting the\n
   *   evaluation of the ith point
   * - The quadrature rule that can be used to initialize the functor on the given element
   */
  using ElemQpArg = std::tuple<const libMesh::Elem *, unsigned int, const QBase *>;

  /**
   * Argument for requesting functor evaluation at quadrature point locations on an element side.
   * Data in the argument:
   * - The element
   * - The element side on which the quadrature points are located
   * - The quadrature point index, e.g. if there are \p n quadrature points, we are requesting the\n
   *   evaluation of the ith point
   * - The quadrature rule that can be used to initialize the functor on the given element and side
   */
  using ElemSideQpArg =
      std::tuple<const libMesh::Elem *, unsigned int, unsigned int, const QBase *>;

  using FunctorType = FunctorImpl<T>;
  using FunctorReturnType = T;
  using ValueType = T;
  using GradientType = typename libMesh::TensorTools::IncrementRank<T>::type;
  using DotType = ValueType;

  virtual ~FunctorImpl() = default;
  FunctorImpl() : _clearance_schedule({EXEC_ALWAYS}) {}

  ///@{
  /**
   * Same as their \p evaluate overloads with the same arguments but allows for caching
   * implementation. These are the methods a user will call in their code
   */
  ValueType operator()(const libMesh::Elem * const & elem, unsigned int state = 0) const;
  ValueType operator()(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const;
  ValueType operator()(const FaceArg & face, unsigned int state = 0) const;
  ValueType operator()(const SingleSidedFaceArg & face, unsigned int state = 0) const;
  ValueType operator()(const ElemQpArg & qp, unsigned int state = 0) const;
  ValueType operator()(const ElemSideQpArg & qp, unsigned int state = 0) const;
  ///@}

  ///@{
  /**
   * Same as their \p evaluateGradient overloads with the same arguments but allows for caching
   * implementation. These are the methods a user will call in their code
   */
  GradientType gradient(const libMesh::Elem * const & elem, unsigned int state = 0) const;
  GradientType gradient(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const;
  GradientType gradient(const FaceArg & face, unsigned int state = 0) const;
  GradientType gradient(const SingleSidedFaceArg & face, unsigned int state = 0) const;
  GradientType gradient(const ElemQpArg & qp, unsigned int state = 0) const;
  GradientType gradient(const ElemSideQpArg & qp, unsigned int state = 0) const;
  ///@}

  ///@{
  /**
   * Same as their \p evaluateDot overloads with the same arguments but allows for caching
   * implementation. These are the methods a user will call in their code
   */
  DotType dot(const libMesh::Elem * const & elem, unsigned int state = 0) const;
  DotType dot(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const;
  DotType dot(const FaceArg & face, unsigned int state = 0) const;
  DotType dot(const SingleSidedFaceArg & face, unsigned int state = 0) const;
  DotType dot(const ElemQpArg & qp, unsigned int state = 0) const;
  DotType dot(const ElemSideQpArg & qp, unsigned int state = 0) const;
  ///@}

  void residualSetup() override;
  void jacobianSetup() override;
  void timestepSetup() override;

  /**
   * Set how often to clear the functor evaluation cache
   */
  void setCacheClearanceSchedule(const std::set<ExecFlagType> & clearance_schedule);

  /**
   * Returns whether this is an extrapolated boundary face
   */
  virtual bool isExtrapolatedBoundaryFace(const FaceInfo &) const { return false; }

protected:
  /**
   * Evaluate the functor with a given element. Some example implementations of this method
   * could compute an element-average or evaluate at the element centroid
   */
  virtual ValueType evaluate(const libMesh::Elem * const & elem, unsigned int state) const = 0;

  /**
   * @param elem_from_face See the \p ElemFromFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const ElemFromFaceArg & elem_from_face, unsigned int state) const = 0;

  /**
   * @param face See the \p FaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const FaceArg & face, unsigned int state) const = 0;

  /**
   * @param face See the \p SingleSidedFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const SingleSidedFaceArg & face, unsigned int state) const = 0;

  /**
   * @param qp See the \p ElemQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const ElemQpArg & qp, unsigned int state) const = 0;

  /**
   * @param side_qp See the \p ElemSideQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor evaluated at the requested time and space
   */
  virtual ValueType evaluate(const ElemSideQpArg & side_qp, unsigned int state) const = 0;

  /**
   * Evaluate the functor gradient with a given element. Some example implementations of this method
   * could compute an element-average or evaluate at the element centroid
   */
  virtual GradientType evaluateGradient(const libMesh::Elem * const &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param elem_from_face See the \p ElemFromFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const ElemFromFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p FaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const FaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p SingleSidedFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const SingleSidedFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param qp See the \p ElemQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const ElemQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param side_qp See the \p ElemSideQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor gradient evaluated at the requested time and space
   */
  virtual GradientType evaluateGradient(const ElemSideQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * Evaluate the functor time derivative with a given element. Some example implementations of this
   * method could compute an element-average or evaluate at the element centroid
   */
  virtual DotType evaluateDot(const libMesh::Elem * const &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param elem_from_face See the \p ElemFromFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const ElemFromFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p FaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const FaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param face See the \p SingleSidedFaceArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const SingleSidedFaceArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param qp See the \p ElemQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const ElemQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

  /**
   * @param side_qp See the \p ElemSideQpArg doxygen
   * @param state Corresponds to a time argument. A value of 0 corresponds to current time, 1
   * corresponds to the old time, 2 corresponds to the older time, etc.
   * @return The functor time derivative evaluated at the requested time and space
   */
  virtual DotType evaluateDot(const ElemSideQpArg &, unsigned int) const
  {
    mooseError("not implemented");
  }

private:
  /**
   * clear cache data
   */
  void clearCacheData();

  /**
   * check a qp cache and if invalid then evaluate
   */
  template <typename SpaceArg, typename TimeArg>
  ValueType queryQpCache(unsigned int qp,
                         const QBase & qrule,
                         std::vector<std::pair<bool, T>> & qp_cache_data,
                         const SpaceArg & space,
                         const TimeArg & time) const;

  /// How often to clear the material property cache
  std::set<ExecFlagType> _clearance_schedule;

  // Data for traditional element-quadrature point property evaluations which are useful for caching
  // implementation

  /// Current key for qp map cache
  mutable dof_id_type _current_qp_map_key = DofObject::invalid_id;

  /// Current value for qp map cache
  mutable std::vector<std::pair<bool, T>> * _current_qp_map_value = nullptr;

  /// Cached element quadrature point functor property evaluations. The map key is the element
  /// id. The map values should have size corresponding to the number of quadrature points on the
  /// element. The vector elements are pairs. The first member of the pair indicates whether a
  /// cached value has been computed. The second member of the pair is the (cached) value. If the
  /// boolean is false, then the value cannot be trusted
  mutable std::unordered_map<dof_id_type, std::vector<std::pair<bool, T>>> _qp_to_value;

  // Data for traditional element-side-quadrature point property evaluations which are useful for
  // caching implementation

  /// Current key for side-qp map cache
  mutable dof_id_type _current_side_qp_map_key = DofObject::invalid_id;

  /// Current value for side-qp map cache
  mutable std::vector<std::vector<std::pair<bool, T>>> * _current_side_qp_map_value = nullptr;

  /// Cached element quadrature point functor property evaluations. The map key is the element
  /// id. The map values are a multi-dimensional vector (or vector of vectors) with the first index
  /// corresponding to the side and the second index corresponding to the quadrature point
  /// index. The elements returned after double indexing are pairs. The first member of the pair
  /// indicates whether a cached value has been computed. The second member of the pair is the
  /// (cached) value. If the boolean is false, then the value cannot be trusted
  mutable std::unordered_map<dof_id_type, std::vector<std::vector<std::pair<bool, T>>>>
      _side_qp_to_value;

  template <typename>
  friend class Functor;
};

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const Elem * const & elem, const unsigned int state) const
{
  return evaluate(elem, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemFromFaceArg & elem_from_face, const unsigned int state) const
{
  return evaluate(elem_from_face, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const FaceArg & face, const unsigned int state) const
{
  return evaluate(face, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const SingleSidedFaceArg & face, const unsigned int state) const
{
  return evaluate(face, state);
}

template <typename T>
template <typename SpaceArg, typename TimeArg>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::queryQpCache(const unsigned int qp,
                             const QBase & qrule,
                             std::vector<std::pair<bool, T>> & qp_cache_data,
                             const SpaceArg & space,
                             const TimeArg & time) const
{
  // Check and see whether we even have sized for this quadrature point. If we haven't then we must
  // evaluate
  if (qp >= qp_cache_data.size())
  {
    qp_cache_data.resize(qrule.n_points(), std::make_pair(false, T()));
    auto & pr = qp_cache_data[qp];
    pr.second = evaluate(space, time);
    pr.first = true;
    return pr.second;
  }

  // We've already sized for this qp, so let's see whether we have a valid cache value
  auto & pr = qp_cache_data[qp];
  if (pr.first)
    return pr.second;

  // No valid cache value so evaluate
  pr.second = evaluate(space, time);
  pr.first = true;
  return pr.second;
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemQpArg & elem_qp, const unsigned int state) const
{
  if (_clearance_schedule.count(EXEC_ALWAYS))
    return evaluate(elem_qp, state);

  const auto elem_id = std::get<0>(elem_qp)->id();
  if (elem_id != _current_qp_map_key)
  {
    _current_qp_map_key = elem_id;
    _current_qp_map_value = &_qp_to_value[elem_id];
  }
  auto & qp_data = *_current_qp_map_value;
  const auto qp = std::get<1>(elem_qp);
  const auto * const qrule = std::get<2>(elem_qp);
  mooseAssert(qrule, "qrule must be non-null");

  return queryQpCache(qp, *qrule, qp_data, elem_qp, state);
}

template <typename T>
typename FunctorImpl<T>::ValueType
FunctorImpl<T>::operator()(const ElemSideQpArg & elem_side_qp, const unsigned int state) const
{
  if (_clearance_schedule.count(EXEC_ALWAYS))
    return evaluate(elem_side_qp, state);

  const Elem * const elem = std::get<0>(elem_side_qp);
  mooseAssert(elem, "elem must be non-null");
  const auto elem_id = elem->id();
  if (elem_id != _current_side_qp_map_key)
  {
    _current_side_qp_map_key = elem_id;
    _current_side_qp_map_value = &_side_qp_to_value[elem_id];
  }
  auto & side_qp_data = *_current_side_qp_map_value;
  const auto side = std::get<1>(elem_side_qp);
  const auto qp = std::get<2>(elem_side_qp);
  const auto * const qrule = std::get<3>(elem_side_qp);
  mooseAssert(qrule, "qrule must be non-null");

  // Check and see whether we even have sized for this side
  if (side >= side_qp_data.size())
    side_qp_data.resize(elem->n_sides());

  // Ok we were sized enough for our side
  auto & qp_data = side_qp_data[side];
  return queryQpCache(qp, *qrule, qp_data, elem_side_qp, state);
}

template <typename T>
void
FunctorImpl<T>::setCacheClearanceSchedule(const std::set<ExecFlagType> & clearance_schedule)
{
  _clearance_schedule = clearance_schedule;
}

template <typename T>
void
FunctorImpl<T>::clearCacheData()
{
  for (auto & map_pr : _qp_to_value)
    for (auto & pr : map_pr.second)
      pr.first = false;

  for (auto & map_pr : _side_qp_to_value)
  {
    auto & side_vector = map_pr.second;
    for (auto & qp_vector : side_vector)
      for (auto & pr : qp_vector)
        pr.first = false;
  }

  _current_qp_map_key = DofObject::invalid_id;
  _current_qp_map_value = nullptr;
  _current_side_qp_map_key = DofObject::invalid_id;
  _current_side_qp_map_value = nullptr;
}

template <typename T>
void
FunctorImpl<T>::timestepSetup()
{
  if (_clearance_schedule.count(EXEC_TIMESTEP_BEGIN))
    clearCacheData();
}

template <typename T>
void
FunctorImpl<T>::residualSetup()
{
  if (_clearance_schedule.count(EXEC_LINEAR))
    clearCacheData();
}

template <typename T>
void
FunctorImpl<T>::jacobianSetup()
{
  if (_clearance_schedule.count(EXEC_NONLINEAR))
    clearCacheData();
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const Elem * const & elem, const unsigned int state) const
{
  return evaluateGradient(elem, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemFromFaceArg & elem_from_face, const unsigned int state) const
{
  return evaluateGradient(elem_from_face, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const FaceArg & face, const unsigned int state) const
{
  return evaluateGradient(face, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const SingleSidedFaceArg & face, const unsigned int state) const
{
  return evaluateGradient(face, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemQpArg & elem_qp, const unsigned int state) const
{
  return evaluateGradient(elem_qp, state);
}

template <typename T>
typename FunctorImpl<T>::GradientType
FunctorImpl<T>::gradient(const ElemSideQpArg & elem_side_qp, const unsigned int state) const
{
  return evaluateGradient(elem_side_qp, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const Elem * const & elem, const unsigned int state) const
{
  return evaluateDot(elem, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemFromFaceArg & elem_from_face, const unsigned int state) const
{
  return evaluateDot(elem_from_face, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const FaceArg & face, const unsigned int state) const
{
  return evaluateDot(face, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const SingleSidedFaceArg & face, const unsigned int state) const
{
  return evaluateDot(face, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemQpArg & elem_qp, const unsigned int state) const
{
  return evaluateDot(elem_qp, state);
}

template <typename T>
typename FunctorImpl<T>::DotType
FunctorImpl<T>::dot(const ElemSideQpArg & elem_side_qp, const unsigned int state) const
{
  return evaluateDot(elem_side_qp, state);
}

/**
 * This is a wrapper that forwards calls to the implementation,
 * which can be switched out at any time without disturbing references to
 * FunctorImpl. Implementation motivated by https://stackoverflow.com/a/65455485/4493669
 */
template <typename T>
class Functor : public FunctorImpl<T>
{
public:
  using typename Moose::FunctorImpl<T>::FaceArg;
  using typename Moose::FunctorImpl<T>::SingleSidedFaceArg;
  using typename Moose::FunctorImpl<T>::ElemFromFaceArg;
  using typename Moose::FunctorImpl<T>::ElemQpArg;
  using typename Moose::FunctorImpl<T>::ElemSideQpArg;
  using typename Moose::FunctorImpl<T>::FunctorType;
  using typename Moose::FunctorImpl<T>::FunctorReturnType;
  using typename Moose::FunctorImpl<T>::ValueType;
  using typename Moose::FunctorImpl<T>::GradientType;
  using typename Moose::FunctorImpl<T>::DotType;

  /**
   * Construct wrapper from wrapped object
   */
  Functor(const FunctorImpl<T> & wrapped) : FunctorImpl<T>(), _wrapped(&wrapped) {}

  template <typename FunctorIn>
  Functor(FunctorIn && wrapped)
    : FunctorImpl<T>(),
      _owned(std::make_unique<FunctorIn>(std::move(wrapped))),
      _wrapped(_owned.get())
  {
  }

  template <typename FunctorIn>
  Functor(std::unique_ptr<FunctorIn> && wrapped)
    : FunctorImpl<T>(), _owned(std::move(wrapped)), _wrapped(_owned.get())
  {
  }

  /**
   * Assign our wrapped object to be something new and release our previously wrapped object
   */
  Functor<T> & operator=(const FunctorImpl<T> & wrapped)
  {
    _owned.reset();
    _wrapped = &wrapped;
    return *this;
  }

  template <typename FunctorIn>
  Functor<T> & operator=(FunctorIn && wrapped)
  {
    _owned.reset();
    _owned = std::make_unique<FunctorIn>(std::move(wrapped));
    _wrapped = _owned.get();
    return *this;
  }

  template <typename FunctorIn>
  Functor<T> & operator=(std::unique_ptr<FunctorIn> && wrapped)
  {
    _owned.reset();
    _owned = std::move(wrapped);
    _wrapped = _owned.get();
    return *this;
  }

  virtual ~Functor() = default;

  /**
   * Tests whether the wrapped object is of the requested type
   */
  template <typename T2>
  bool wrapsType() const
  {
    return dynamic_cast<const T2 *>(_wrapped);
  }

protected:
  ///@{
  /**
   * Forward evaluate calls to wrapped object
   */
  ValueType evaluate(const libMesh::Elem * const & elem, unsigned int state = 0) const override
  {
    return _wrapped->evaluate(elem, state);
  }
  ValueType evaluate(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const override
  {
    return _wrapped->evaluate(elem_from_face, state);
  }
  ValueType evaluate(const FaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->evaluate(face, state);
  }
  ValueType evaluate(const SingleSidedFaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->evaluate(face, state);
  }
  ValueType evaluate(const ElemQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->evaluate(qp, state);
  }
  ValueType evaluate(const ElemSideQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->evaluate(qp, state);
  }

  GradientType evaluateGradient(const libMesh::Elem * const & elem,
                                unsigned int state = 0) const override
  {
    return _wrapped->evaluateGradient(elem, state);
  }
  GradientType evaluateGradient(const ElemFromFaceArg & elem_from_face,
                                unsigned int state = 0) const override
  {
    return _wrapped->evaluateGradient(elem_from_face, state);
  }
  GradientType evaluateGradient(const FaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->evaluateGradient(face, state);
  }
  GradientType evaluateGradient(const SingleSidedFaceArg & face,
                                unsigned int state = 0) const override
  {
    return _wrapped->evaluateGradient(face, state);
  }
  GradientType evaluateGradient(const ElemQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->evaluateGradient(qp, state);
  }
  GradientType evaluateGradient(const ElemSideQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->evaluateGradient(qp, state);
  }

  DotType evaluateDot(const libMesh::Elem * const & elem, unsigned int state = 0) const override
  {
    return _wrapped->evaluateDot(elem, state);
  }
  DotType evaluateDot(const ElemFromFaceArg & elem_from_face, unsigned int state = 0) const override
  {
    return _wrapped->evaluateDot(elem_from_face, state);
  }
  DotType evaluateDot(const FaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->evaluateDot(face, state);
  }
  DotType evaluateDot(const SingleSidedFaceArg & face, unsigned int state = 0) const override
  {
    return _wrapped->evaluateDot(face, state);
  }
  DotType evaluateDot(const ElemQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->evaluateDot(qp, state);
  }
  DotType evaluateDot(const ElemSideQpArg & qp, unsigned int state = 0) const override
  {
    return _wrapped->evaluateDot(qp, state);
  }
  ///@}

  /// Our wrapped object
  std::unique_ptr<FunctorImpl<T>> _owned;
  const FunctorImpl<T> * _wrapped;
};

/**
 * Class template for creating constants
 */
template <typename T>
class ConstantFunctor : public FunctorImpl<T>
{
public:
  using typename FunctorImpl<T>::FaceArg;
  using typename FunctorImpl<T>::SingleSidedFaceArg;
  using typename FunctorImpl<T>::ElemFromFaceArg;
  using typename FunctorImpl<T>::ElemQpArg;
  using typename FunctorImpl<T>::ElemSideQpArg;
  using typename FunctorImpl<T>::FunctorType;
  using typename FunctorImpl<T>::FunctorReturnType;
  using typename FunctorImpl<T>::ValueType;
  using typename FunctorImpl<T>::GradientType;
  using typename FunctorImpl<T>::DotType;

  ConstantFunctor(const ValueType & value) : _value(value) {}
  ConstantFunctor(ValueType && value) : _value(value) {}

private:
  ValueType evaluate(const libMesh::Elem * const &, unsigned int) const override final
  {
    return _value;
  }
  ValueType evaluate(const ElemFromFaceArg &, unsigned int) const override final { return _value; }
  ValueType evaluate(const FaceArg &, unsigned int) const override final { return _value; }
  ValueType evaluate(const SingleSidedFaceArg &, unsigned int) const override final
  {
    return _value;
  }
  ValueType evaluate(const ElemQpArg &, unsigned int) const override final { return _value; }
  ValueType evaluate(const ElemSideQpArg &, unsigned int) const override final { return _value; }

  GradientType evaluateGradient(const libMesh::Elem * const &, unsigned int) const override final
  {
    return 0;
  }
  GradientType evaluateGradient(const ElemFromFaceArg &, unsigned int) const override final
  {
    return 0;
  }
  GradientType evaluateGradient(const FaceArg &, unsigned int) const override final { return 0; }
  GradientType evaluateGradient(const SingleSidedFaceArg &, unsigned int) const override final
  {
    return 0;
  }
  GradientType evaluateGradient(const ElemQpArg &, unsigned int) const override final { return 0; }
  GradientType evaluateGradient(const ElemSideQpArg &, unsigned int) const override final
  {
    return 0;
  }

  DotType evaluateDot(const libMesh::Elem * const &, unsigned int) const override final
  {
    return 0;
  }
  DotType evaluateDot(const ElemFromFaceArg &, unsigned int) const override final { return 0; }
  DotType evaluateDot(const FaceArg &, unsigned int) const override final { return 0; }
  DotType evaluateDot(const SingleSidedFaceArg &, unsigned int) const override final { return 0; }
  DotType evaluateDot(const ElemQpArg &, unsigned int) const override final { return 0; }
  DotType evaluateDot(const ElemSideQpArg &, unsigned int) const override final { return 0; }

private:
  ValueType _value;
};
}
