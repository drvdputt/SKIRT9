/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SELFCONSISTENTOPACITYOPTIONS_HPP
#define SELFCONSISTENTOPACITYOPTIONS_HPP

#include "ItemInfo.hpp"
#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** The SelfConsistentOpacityOptions class simply offers a number of configuration options related
    to the self-consistent calculation of the dust and gas opacity, including convergence criteria
    for the iteration process. These options are relevant only when self-consistent opacity
    iterations are enabled for the simulations, and a medium with a radiation field-dependent
    opacity is present. */
class SelfConsistentOpacityOptions : public SimulationItem
{
    ITEM_CONCRETE(SelfConsistentOpacityOptions, SimulationItem,
                  "a set of options related to self-consistently calculating the opacity")

        PROPERTY_BOOL(withSecondary,
                      "take secondary emission into account during the self-consistent opacity iteration")
        ATTRIBUTE_DEFAULT_VALUE(withSecondary, "true")
        ATTRIBUTE_DISPLAYED_IF(withSecondary, "Level3")

        PROPERTY_INT(minIterations, "the minimum number of opacity iterations")
        ATTRIBUTE_MIN_VALUE(minIterations, "2")
        ATTRIBUTE_MAX_VALUE(minIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(minIterations, "3")
        ATTRIBUTE_DISPLAYED_IF(minIterations, "Level3")

        PROPERTY_INT(maxIterations, "the maximum number of opacity iterations")
        ATTRIBUTE_MIN_VALUE(maxIterations, "2")
        ATTRIBUTE_MAX_VALUE(maxIterations, "1000")
        ATTRIBUTE_DEFAULT_VALUE(maxIterations, "10")
        ATTRIBUTE_DISPLAYED_IF(maxIterations, "Level3")

        PROPERTY_DOUBLE(maxFractionOfPrimary_gas,
                        "gas convergence is reached when the total gas-absorbed secondary luminosity  "
                        "is less than this fraction of the gas absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrimary_gas, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrimary_gas, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrimary_gas, "0.01")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrimary_gas, "Level2")

        PROPERTY_DOUBLE(maxFractionOfPrimary_dust,
                        "dust convergence is reached when the total dust-absorbed secondary luminosity  "
                        "is less than this fraction of the dust absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrimary_dust, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrimary_dust, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrimary_dust, "0.01")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrimary_dust, "Level2")

        PROPERTY_DOUBLE(
            maxFractionOfPrevious_gas,
            "gas convergence is reached when the total gas-absorbed luminosity has changed by less than this "
            "fraction compared to the previous iteration for both the primary and the secondary emission")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrevious_gas, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrevious_gas, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrevious_gas, "0.03")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrevious_gas, "Level2")

        PROPERTY_DOUBLE(
            maxFractionOfPrevious_dust,
            "dust convergence is reached when the total dust-absorbed luminosity has changed by less than this "
            "fraction compared to the previous iteration for both the primary and the secondary emission")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPrevious_dust, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPrevious_dust, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPrevious_dust, "0.03")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPrevious_dust, "Level2")

        PROPERTY_DOUBLE(primaryPacketsMultiplier,
                        "the multiplier on the number of primary photon packets lauched for each opacity iteration")
        ATTRIBUTE_MIN_VALUE(primaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(primaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(primaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(primaryPacketsMultiplier, "Level3")

        PROPERTY_DOUBLE(secondaryPacketsMultiplier,
                        "the multiplier on the number of secondary photon packets launched for each opacity iteration")
        ATTRIBUTE_MIN_VALUE(secondaryPacketsMultiplier, "]0")
        ATTRIBUTE_MAX_VALUE(secondaryPacketsMultiplier, "1000]")
        ATTRIBUTE_DEFAULT_VALUE(secondaryPacketsMultiplier, "1")
        ATTRIBUTE_DISPLAYED_IF(secondaryPacketsMultiplier, "Level3")

    ITEM_END()
};

////////////////////////////////////////////////////////////////////

#endif
