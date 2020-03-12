/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GasProbe.hpp"
#include "Gas.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "TextOutFile.hpp"

////////////////////////////////////////////////////////////////////

void GasProbe::probeRun()
{
    auto ms = find<MediumSystem>(false);
    auto grid = ms->grid();
    auto log = find<Log>(false);
    if (!ms->hasGas())
    {
        log->warning("No gas is present! Gas probe will not run!");
        return;
    }
    TextOutFile file(this, "gas", "gas properties per cell");
    for (auto s : {"index", "x", "y", "z", "T", "np", "nH", "nH2"}) file.addColumn(s, "", 'd');
    int numCells = grid->numCells();
    for (int m = 0; m < numCells; m++)
    {
        Position p = grid->centralPositionInCell(m);
        file.writeRow(vector<double>{static_cast<double>(m), p.x(), p.y(), p.z(), Gas::temperature(m), Gas::np(m),
                                     Gas::nH(m), Gas::nH2(m)});
    }
}
