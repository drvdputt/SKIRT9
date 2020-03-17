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
    for (auto& s : {"x", "y", "z"}) file.addColumn(s, "m", 'e');
    int numCols = 4;

    // temperature
    file.addColumn("T", "K", 'e');
    numCols++;

    // abundances
    for (auto& s : {"np", "nH", "nH2"}) file.addColumn(s, "cm-3", 'e');
    numCols += 3;

    if (extendedDiagnostics())
    {
        // TODO: units
        for (auto& s : Gas::diagnosticNames())
        {
            file.addColumn(s, "", 'e');
            numCols++;
        }
    }

    // do the diagnostics calculation in parallel (mostly relevent for the extended option)
    int numCells = grid->numCells();
    Table<2> numbers(numCells, numCols);
    find<ParallelFactory>()->parallelDistributed()->call(numCells, [&](size_t firstIndex, size_t numIndices) {
        for (size_t m = firstIndex; m < firstIndex + numIndices; m++)
        {
            vector<double> numbersForCell;
            numbersForCell.reserve(numCols);

            // basic properties
            Position p = grid->centralPositionInCell(m);
            for (auto value : {static_cast<double>(m), p.x(), p.y(), p.z(), Gas::temperature(m), Gas::np(m), Gas::nH(m),
                               Gas::nH2(m)})
                numbersForCell.push_back(value);

            if (extendedDiagnostics())
            {
                double n = 0;
                for (int h = 0; h != ms->numMedia(); ++h)
                    if (ms->isGas(h)) n += ms->numberDensity(m, h);

                auto hv = Gas::hIndices();
                Array nv(hv.size());
                int c = 0;
                for (int h : hv)
                {
                    nv[c] = ms->numberDensity(m, h);
                    c++;
                }

                auto diagnostics = Gas::diagnostics(m, n, ms->meanIntensity(m), nv);
                for (auto value : diagnostics) numbersForCell.push_back(value);
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
