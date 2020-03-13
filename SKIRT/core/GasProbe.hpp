/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GASPROBE_HPP
#define GASPROBE_HPP

#include "ItemInfo.hpp"
#include "Probe.hpp"

////////////////////////////////////////////////////////////////////

/** GasProbe outputs a column text files (named <tt>prefix_probe_gas.dat</tt>) listing the gas
    temperature and abundances throughout the model grid. The output file contains a line for each
    cell in the spatial grid. Each line contains the cell index and the coordinates of the center
    of the cell, followed by the properties. When \c extendedDiagnostics is enabled, the
    equilibrium calculation will be re-done, but with the diagnostic option enabled (see \c
    GasDiagnostics* argument of \c GasModule::GasInterface::updateGasState). Many more columns will
    be added to the output file, including (but not limited to) the heating and cooling
    contributions, formation and destruction rates of certain species, and the adjusted grain
    temperature used internally by the gas code. */
class GasProbe : public Probe
{
    ITEM_CONCRETE(GasProbe, Probe, "information about the gas at the end of the simulation")
        ATTRIBUTE_TYPE_DISPLAYED_IF(GasProbe, "GasMedium")

        PROPERTY_BOOL(extendedDiagnostics, "enable advanced gas diagnostics (potentially slow)")
        ATTRIBUTE_DEFAULT_VALUE(extendedDiagnostics, "false")

    ITEM_END()

    //======================== Other Functions =======================

public:
    /** This function performs probing after all photon packets have been emitted and detected */
    void probeRun() override;
};

////////////////////////////////////////////////////////////////////

#endif
