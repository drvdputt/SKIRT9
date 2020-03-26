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

        PROPERTY_DOUBLE(FractionOfPrimaryGas,
                        "gas convergence is reached when the total gas-absorbed secondary luminosity "
                        "is less than this fraction of the gas-absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(FractionOfPrimaryGas, "]0")
        ATTRIBUTE_MAX_VALUE(FractionOfPrimaryGas, "1[")
        ATTRIBUTE_DEFAULT_VALUE(FractionOfPrimaryGas, "0.01")
        ATTRIBUTE_DISPLAYED_IF(FractionOfPrimaryGas, "Level2")

        PROPERTY_DOUBLE(FractionOfPrimaryDust,
                        "dust convergence is reached when the total dust-absorbed secondary luminosity "
                        "is less than this fraction of the dust-absorbed primary luminosity")
        ATTRIBUTE_MIN_VALUE(FractionOfPrimaryDust, "]0")
        ATTRIBUTE_MAX_VALUE(FractionOfPrimaryDust, "1[")
        ATTRIBUTE_DEFAULT_VALUE(FractionOfPrimaryDust, "0.01")
        ATTRIBUTE_DISPLAYED_IF(FractionOfPrimaryDust, "Level2")

        PROPERTY_DOUBLE(
            maxFractionOfPreviousGas,
            "gas convergence is reached when the total gas-absorbed luminosity has changed by less than this "
            "fraction compared to the previous iteration for both the primary and the secondary emission")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPreviousGas, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPreviousGas, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPreviousGas, "0.03")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPreviousGas, "Level2")

        PROPERTY_DOUBLE(
            maxFractionOfPreviousDust,
            "dust convergence is reached when the total dust-absorbed luminosity has changed by less than this "
            "fraction compared to the previous iteration for both the primary and the secondary emission")
        ATTRIBUTE_MIN_VALUE(maxFractionOfPreviousDust, "]0")
        ATTRIBUTE_MAX_VALUE(maxFractionOfPreviousDust, "1[")
        ATTRIBUTE_DEFAULT_VALUE(maxFractionOfPreviousDust, "0.03")
        ATTRIBUTE_DISPLAYED_IF(maxFractionOfPreviousDust, "Level2")

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
