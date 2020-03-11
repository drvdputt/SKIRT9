/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASPROBE_HPP
#define GASPROBE_HPP

#include "ItemInfo.hpp"
#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** GasProbe outputs a column text files (named <tt>prefix_probe_gas.dat</tt>) listing temperature
    and abundances (and more) gas properties throughout the model. The output file contains a line
    for each cell in the spatial grid. Each line contains the cell index and the coordinates of the
    center of the cell, followed by the properties. */
class GasProbe : public Probe
{
    ITEM_CONCRETE(GasProbe, Probe, "information about the gas at the end of the simulation")
        ATTRIBUTE_TYPE_DISPLAYED_IF(GasProbe, "GasMedium")

    ITEM_END()

    //======================== Other Functions =======================

    public:

    /** This function performs probing after all photon packets have been emitted and detected */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
