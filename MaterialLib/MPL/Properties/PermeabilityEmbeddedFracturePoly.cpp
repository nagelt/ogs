/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/PermeabilityEmbeddedFracturePoly.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/KelvinVector.h"

#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
template <int DisplacementDim>
PermeabilityEmbeddedFracturePoly<DisplacementDim>::PermeabilityEmbeddedFracturePoly(
    double const intrinsic_permeability,
    std::vector<double> const mean_fracture_distances,
    std::vector<double> const threshold_strains,
    Eigen::Matrix<double, 3, 3> const fracture_normals,
    ParameterLib::Parameter<double> const& fracture_rotation)
    : _n(fracture_normals),
      _k(intrinsic_permeability),
      _a(mean_fracture_distances),
      _e0(threshold_strains),
      _phi(fracture_rotation),
      _b0(std::sqrt(12 * _k))
{};

template <int DisplacementDim>
PropertyDataType PermeabilityEmbeddedFracturePoly<DisplacementDim>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    constexpr int symmetric_tensor_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using SymmetricTensor = Eigen::Matrix<double, symmetric_tensor_size, 1>;

    auto const eps = formEigenTensor<3>(std::get<SymmetricTensor>(
        variable_array[static_cast<int>(Variable::strain)]));

    SymmetricTensor result = SymmetricTensor::Zero();

    double const phi = std::get<double>(fromVector(_phi(t, pos)));

    Eigen::Matrix3d const rotMat = (Eigen::Matrix3d() << cos(phi), -sin(phi), 0,
                                    sin(phi), cos(phi), 0, 0, 0, 1)
                                       .finished();

    for (int i = 0; i < 3; i++)
    {
        Eigen::Matrix<double, 3, 1> const ni = rotMat * _n.col(i);
        double const e_n = (eps * ni).dot(ni.transpose());

        double const H_de = (e_n > _e0[i]) ? 1.0 : 0.0;
        double const b_f = _b0 + H_de * _a[i] * (e_n - _e0[i]);
        double const coeff = H_de * (b_f / _a[i]) * ((b_f * b_f / 12.0) - _k);

        SymmetricTensor part_result;
        if (DisplacementDim == 2)
        {
            part_result << 1 - ni[0] * ni[0], 1 - ni[1] * ni[1], 1 - ni[2] * ni[2],
                ni[0] * ni[1];
        }
        else if (DisplacementDim == 3)
        {
            part_result << 1 - ni[0] * ni[0], 1 - ni[1] * ni[1], 1 - ni[2] * ni[2],
                ni[0] * ni[1], ni[1] * ni[2], ni[0] * ni[2];
        }
        part_result *= coeff;
        result += part_result;
    }

    SymmetricTensor k_I = SymmetricTensor::Zero();
    k_I.template head<3>() = Eigen::Vector3d::Constant(_k);
    result += k_I;

    // _k * I + coeff * (I - M)
    return result;
}
template <int DisplacementDim>
PropertyDataType PermeabilityEmbeddedFracturePoly<DisplacementDim>::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::strain) &&
           "PermeabilityEmbeddedFracturePoly::dValue is implemented for "
           " derivatives with respect to strain only.");

    constexpr int symmetric_tensor_size =
        MathLib::KelvinVector::KelvinVectorDimensions<DisplacementDim>::value;
    using SymmetricTensor = Eigen::Matrix<double, symmetric_tensor_size, 1>;

    auto const eps = formEigenTensor<3>(std::get<SymmetricTensor>(
        variable_array[static_cast<int>(Variable::strain)]));

    SymmetricTensor result = SymmetricTensor::Zero();
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();

    for (int i = 0; i < 3; i++)
    {
        Eigen::Matrix<double, 3, 1> const ni = _n.col(i);
        double const e_n = (eps * ni).dot(ni.transpose());

        double const H_de = (e_n > _e0[i]) ? 1.0 : 0.0;
        double const b_f = _b0 + H_de * _a[i] * (e_n - _e0[i]);

        Eigen::Matrix3d const M = ni * ni.transpose();
        Eigen::Matrix3d const dkde = H_de * (b_f * b_f / 4 - _k) * (I - M) * M;

        SymmetricTensor part_result;

        if (DisplacementDim == 2)
        {
            part_result << dkde(0, 0), dkde(1, 1), dkde(2, 2), dkde(0, 1);
        }
        else if (DisplacementDim == 3)
        {
            part_result << dkde(0, 0), dkde(1, 1), dkde(2, 2), dkde(0, 1),
                dkde(1, 2), dkde(0, 2);
        }
        result += part_result;
    }

    //auto const dk_de = H_de * (b_f * b_f / 4 - k_m) * (I - M) * M;
    return result;
}

template class PermeabilityEmbeddedFracturePoly<2>;
template class PermeabilityEmbeddedFracturePoly<3>;
}  // namespace MaterialPropertyLib
