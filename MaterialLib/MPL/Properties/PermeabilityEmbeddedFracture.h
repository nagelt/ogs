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

namespace MaterialPropertyLib
{
class Medium;
class Phase;
class Component;
/**
 * \class PermeabilityEmbeddedFracture
 * \brief Permeability function proposed by Olivella&Alonso
 * \details This property must be a medium property, it
 * computes the permeability in dependence of the strain
 */
template <int DisplacementDim>
class PermeabilityEmbeddedFracture final : public Property
{
private:
    Medium* _medium = nullptr;
    double const _k;
    double const _a;
    double const _e0;
    double const _b0;
    bool const _n_const;
    Eigen::Matrix<double, 3, 1> const _n;
 public:
    PermeabilityEmbeddedFracture(Eigen::Matrix<double, 3, 1> const fracture_normal,
                         bool const fracture_normal_is_constant,
                         double const intrinsic_permeability,
                         double const mean_fracture_distance,
                         double const threshold_strain);
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
                "The property 'PermeabilityEmbeddedFracture' is implemented on the "
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
