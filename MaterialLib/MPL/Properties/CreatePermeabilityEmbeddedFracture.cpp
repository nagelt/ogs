/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "PermeabilityEmbeddedFracture.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createPermeabilityEmbeddedFracture(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "PermeabilityEmbeddedFracture");
    DBUG("Create PermeabilityEmbeddedFracture medium property");

    double const intrinsic_permeability =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracture__intrinsic_permeability}
        config.getConfigParameter<double>("intrinsic_permeability");

    double const mean_fracture_distance =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracture__mean_frac_distance}
        config.getConfigParameter<double>("mean_frac_distance");

    double const threshold_strain =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracture__threshold_strain}
        config.getConfigParameter<double>("threshold_strain");

    bool fracture_normal_is_constant = false;
    Eigen::Matrix<double, 3, 1> fracture_normal;
    //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracture__fracture_normal}
    if (auto const n_ptr =
            config.getConfigParameterOptional<std::vector<double>>(
                "fracture_normal"))
    {
        if ((*n_ptr).size() != 3)
        {
            OGS_FATAL(
                "The size of the fracture normal vector must be 3, but is %d.",
                (*n_ptr).size());
        }
        DBUG("Using constant fracture normal vector.");
        std::copy_n((*n_ptr).data(), 3, fracture_normal.data());
        fracture_normal_is_constant = true;
    }
    else
    {
        DBUG("No constant fracture normal was given. By default it will be "
             "determined as the third principal stress vector.");
    }

    int const Dim =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracture__Dim}
        config.getConfigParameter<int>("dim");

    if (Dim == 2)
    {
        return std::make_unique<PermeabilityEmbeddedFracture<2>>(
            fracture_normal,
            fracture_normal_is_constant,
            intrinsic_permeability,
            mean_fracture_distance,
            threshold_strain);
    }
    if (Dim == 3)
    {
        return std::make_unique<PermeabilityEmbeddedFracture<3>>(
            fracture_normal,
            fracture_normal_is_constant,
            intrinsic_permeability,
            mean_fracture_distance,
            threshold_strain);
    }

    OGS_FATAL("The Dimension can only be 2 or 3.");
    }
}  // namespace MaterialPropertyLib
