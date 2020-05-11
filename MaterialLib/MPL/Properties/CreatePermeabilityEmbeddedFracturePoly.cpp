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
#include "ParameterLib/Utils.h"
#include "PermeabilityEmbeddedFracturePoly.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createPermeabilityEmbeddedFracturePoly(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "PermeabilityEmbeddedFracturePoly");
    DBUG("Create PermeabilityEmbeddedFracturePoly medium property");

    double const intrinsic_permeability =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracturePoly__intrinsic_permeability}
        config.getConfigParameter<double>("intrinsic_permeability");

    auto const mean_fracture_distances =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracturePoly__mean_frac_distances}
        config.getConfigParameter<std::vector<double>>("mean_frac_distances");
    if (mean_fracture_distances.size() != 3)
    {
        OGS_FATAL(
            "The size of the mean fracture distances vector must be 3, but is "
            "%d.",
            mean_fracture_distances.size());
    }

    auto const threshold_strains =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracturePoly__threshold_strains}
        config.getConfigParameter<std::vector<double>>("threshold_strains");
    if (threshold_strains.size() != 3)
    {
        OGS_FATAL(
            "The size of the mean threshold strains vector must be 3, but is "
            "%d.",
            threshold_strains.size());
    }

    Eigen::Matrix<double, 3, 3> fracture_normals;
    Eigen::Matrix<double, 3, 1> n1, n2, n3;
    auto const n =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracturePoly__fracture_normals}
        config.getConfigParameter<std::vector<double>>("fracture_normals");
    if (n.size() != 6)
    {
        OGS_FATAL(
            "The size of the fracture normals vector must be 6, but is %d.",
            n.size());
    }
    n1 << n[0], n[1], n[2];
    n2 << n[3], n[4], n[5];
    n1 /= n1.norm();
    n2 /= n2.norm();

    if (n1.dot(n2) > std::numeric_limits<double>::epsilon())
    {
        OGS_FATAL(
            "The given fracture normals are not orthogonal. Please provide two "
            "orthogonal fracture normals",
            n.size());
    }

    n3 = n1.cross(n2);
    fracture_normals << n1[0], n1[1], n1[2],
                        n2[0], n2[1], n2[2],
                        n3[0], n3[1], n3[2];

    std::string const parameter_name =
        //! \ogs_file_param{properties__property__TransportPorosityFromMassBalance__fracture_rotation}
        config.getConfigParameter<std::string>("fracture_rotation");

    auto const& fracture_rotation = ParameterLib::findParameter<double>(
            parameter_name, parameters, 0, nullptr);

    int const Dim =
        //! \ogs_file_param{properties__property__PermeabilityEmbeddedFracturePoly__Dim}
        config.getConfigParameter<int>("dim");

    if (Dim == 2)
    {
        return std::make_unique<PermeabilityEmbeddedFracturePoly<2>>(
            intrinsic_permeability,
            mean_fracture_distances,
            threshold_strains,
            fracture_normals,
            fracture_rotation);
    }
    if (Dim == 3)
    {
        return std::make_unique<PermeabilityEmbeddedFracturePoly<3>>(
            intrinsic_permeability,
            mean_fracture_distances,
            threshold_strains,
            fracture_normals,
            fracture_rotation);
    }

    OGS_FATAL("The Dimension can only be 2 or 3.");
}
}  // namespace MaterialPropertyLib
