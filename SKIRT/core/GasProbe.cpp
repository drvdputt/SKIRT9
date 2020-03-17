/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GasProbe.hpp"
#include "FatalError.hpp"
#include "Gas.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "ProcessManager.hpp"
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

    // cell coordinates
    file.addColumn("index", "", 'd');
    for (auto s : {"x", "y", "z"}) file.addColumn(s, "m", 'e');
    int numCols = 4;

    // temperature
    file.addColumn("T", "K", 'e');
    numCols++;

    // abundances
    for (auto s : {"np", "nH", "nH2"}) file.addColumn(s, "cm-3", 'e');
    numCols += 3;

    if (extendedDiagnostics())
    {
        // add more columns
    }

    // do the diagnostics calculation in parallel (mostly relevent for the extended option)
    int numCells = grid->numCells();
    Table<2> numbers(numCells, numCols);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [&](size_t firstIndex, size_t numIndices) {
        for (size_t m = firstIndex; m < firstIndex + numIndices; m++)
        {
            Position p = grid->centralPositionInCell(m);
            vector<double> numbersForCell = {static_cast<double>(m), p.x(),      p.y(),      p.z(),
                                             Gas::temperature(m),    Gas::np(m), Gas::nH(m), Gas::nH2(m)};

            if (extendedDiagnostics())
            {
                // add more numbers to vector
            }

            // temporary check for typo's
            int numNumbers = numbersForCell.size();
            if (numNumbers != numCols) FATALERROR("Incorrect number of elements for row");
            std::copy(std::begin(numbersForCell), std::end(numbersForCell), &numbers(m, 0));
        }
    });

    // write (root only)
    ProcessManager::sumToRoot(numbers.data());
    for (int m = 0; m < numCells; m++) file.writeRow(vector<double>(&numbers(m, 0), &numbers(m, numCols)));
}

////////////////////////////////////////////////////////////////////
