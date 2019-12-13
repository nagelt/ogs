/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ParameterLib/Parameter.h"

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
struct MechanicsBase;
}
}  // namespace MaterialLib
namespace ProcessLib
{
namespace TH2M
{
template <int DisplacementDim>
struct TH2MProcessData
{
    MeshLib::PropertyVector<int> const* const material_ids = nullptr;

    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map = nullptr;

    /// The constitutive relation for the mechanical part.
    /// \note Linear elasticity is the only supported one in the moment.
    std::map<
        int,
        std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
        solid_materials;
    ParameterLib::Parameter<double> const& reference_temperature;
    /// Specific body forces applied to solid and fluid.
    /// It is usually used to apply gravitational forces.
    /// A vector of displacement dimension's length.
    Eigen::Matrix<double, DisplacementDim, 1> const specific_body_force;

    MeshLib::PropertyVector<double>* gas_pressure_interpolated = nullptr;
    MeshLib::PropertyVector<double>* capillary_pressure_interpolated = nullptr;
    MeshLib::PropertyVector<double>* temperature_interpolated = nullptr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

}  // namespace TH2M
}  // namespace ProcessLib
