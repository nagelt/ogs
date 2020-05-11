/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Property.h"
#include "ParameterLib/Parameter.h"

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class PermeabilityEmbeddedFracturePoly
 * \brief Permeability function proposed by Olivella&Alonso
 * \details This property must be a medium property, it
 * computes the permeability in dependence of the strain
 */
template <int DisplacementDim>
class PermeabilityEmbeddedFracturePoly final : public Property
{
private:
    Medium* _medium = nullptr;
    double const _k;
    std::vector<double> const _a;
    std::vector<double> const _e0;
    double const _b0;
    Eigen::Matrix<double, 3, 3> const _n;
    ParameterLib::Parameter<double> const& _phi;
 public:
    PermeabilityEmbeddedFracturePoly(
        double const intrinsic_permeability,
        std::vector<double> const mean_fracture_distances,
        std::vector<double> const threshold_strains,
        Eigen::Matrix<double, 3, 3> const fracture_normals,
         ParameterLib::Parameter<double> const& fracture_rotation);
    /// This method assigns a pointer to the material object that is the owner
    /// of this property
    void setScale(
        std::variant<Medium*, Phase*, Component*> scale_pointer) override
    {
        if (std::holds_alternative<Medium*>(scale_pointer))
        {
            _medium = std::get<Medium*>(scale_pointer);
        }
        else
        {
            OGS_FATAL(
                "The property 'PermeabilityEmbeddedFracturePoly' is implemented on the "
                "'media' scale only.");
        }
    }

    /// Those methods override the base class implementations and
    /// actually compute and set the property _values and _dValues.
    PropertyDataType value(VariableArray const& variable_array,
                           ParameterLib::SpatialPosition const& pos,
                           double const t, double const dt) const override;
    PropertyDataType dValue(VariableArray const& variable_array,
                            Variable const variable,
                            ParameterLib::SpatialPosition const& pos,
                            double const t, double const dt) const override;
};

}  // namespace MaterialPropertyLib
