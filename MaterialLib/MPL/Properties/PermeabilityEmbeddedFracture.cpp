/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/PermeabilityEmbeddedFracture.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/KelvinVector.h"

#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
template <int DisplacementDim>
PermeabilityEmbeddedFracture<DisplacementDim>::PermeabilityEmbeddedFracture(
    Eigen::Matrix<double, 3, 1> const fracture_normal,
    bool const fracture_normal_is_constant,
    double const intrinsic_permeability,
    double const mean_fracture_distance,
    double const threshold_strain)
    : _n(fracture_normal),
      _n_const(fracture_normal_is_constant),
      _k(intrinsic_permeability),
      _a(mean_fracture_distance),
      _e0(threshold_strain),
      _b0(std::sqrt(12 * _k))
{};

template <int DisplacementDim>
PropertyDataType PermeabilityEmbeddedFracture<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    constexpr int symmetric_tensor_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using SymmetricTensor = Eigen::Matrix<double, symmetric_tensor_size, 1>;

    Eigen::Matrix<double, 3, 1> const n = [&] {
        if (_n_const)
        {
            return _n;
        }
        auto const sigma = formEigenTensor<3>(std::get<SymmetricTensor>(
            variable_array[static_cast<int>(Variable::stress)]));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma);
        return (Eigen::Matrix<double, 3, 1>) e_s.eigenvectors().col(2);
    }();

    auto const eps = formEigenTensor<3>(std::get<SymmetricTensor>(
        variable_array[static_cast<int>(Variable::strain)]));
    double const e_n = (eps * n).dot(n.transpose());

    double const H_de = (e_n > _e0) ? 1.0 : 0.0;
    double const b_f = _b0 + H_de * _a * (e_n - _e0);
    double const coeff = H_de * (b_f / _a) * ((b_f * b_f / 12.0) - _k);

    SymmetricTensor result;
    SymmetricTensor k_I = SymmetricTensor::Zero();
    k_I.template head<3>() = Eigen::Vector3d::Constant(_k);

    if (DisplacementDim == 2)
    {
        result << 1 - n[0] * n[0], 1 - n[1] * n[1], 1 - n[2] * n[2],
            -n[0] * n[1];
    }
    else if (DisplacementDim == 3)
    {
        result << 1 - n[0] * n[0], 1 - n[1] * n[1], 1 - n[2] * n[2],
            -n[0] * n[1], -n[1] * n[2], -n[0] * n[2];
    }
    result *= coeff;
    result += k_I;

    // _k * I + coeff * (I - M)
    return result;
}
template <int DisplacementDim>
PropertyDataType PermeabilityEmbeddedFracture<DisplacementDim>::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::strain) &&
           "PermeabilityEmbeddedFracture::dValue is implemented for "
           " derivatives with respect to strain only.");

    constexpr int symmetric_tensor_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using SymmetricTensor = Eigen::Matrix<double, symmetric_tensor_size, 1>;

    Eigen::Matrix<double, 3, 1> const n = [&] {
        if (_n_const)
        {
            return _n;
        }
        auto const sigma = formEigenTensor<3>(std::get<SymmetricTensor>(
            variable_array[static_cast<int>(Variable::stress)]));
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 3, 3>> e_s(sigma);
        return (Eigen::Matrix<double, 3, 1>) e_s.eigenvectors().col(2);
    }();

    auto const eps = formEigenTensor<3>(std::get<SymmetricTensor>(
        variable_array[static_cast<int>(Variable::strain)]));
    double const e_n = (eps * n).dot(n.transpose());

    double const H_de = (e_n > _e0) ? 1.0 : 0.0;
    double const b_f = _b0 + H_de * _a * (e_n - _e0);

    Eigen::Matrix3d const M = n * n.transpose();
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d const dkde = H_de * (b_f * b_f / 4 - _k) * (I - M) * M;

    SymmetricTensor result;

    if (DisplacementDim == 2)
    {
        result << dkde(0, 0), dkde(1, 1), dkde(2, 2), dkde(0, 1);
    }
    else if (DisplacementDim == 3)
    {
        result << dkde(0, 0), dkde(1, 1), dkde(2, 2), dkde(0, 1), dkde(1, 2),
            dkde(0, 2);
    }

    //auto const dk_de = H_de * (b_f * b_f / 4 - k_m) * (I - M) * M;
    return result;
}

template class PermeabilityEmbeddedFracture<2>;
template class PermeabilityEmbeddedFracture<3>;
}  // namespace MaterialPropertyLib
