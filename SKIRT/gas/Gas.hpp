/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */
#ifndef GAS_HPP
#define GAS_HPP

#include "Array.hpp"

/** The functions defined here interface directly with the external gas module code. The main
    reason for introducing this extra layer of abstraction is that SKIRT can be compiled without
    the gas module being present on the system. The gas properties are accessed and updated through
    a set of static functions. */
class Gas
{
public:
    struct DustInfo
    {
        // An index of choice for bookkeeping. Will be returned by \c dustIndices to indicate which
        // densities are needed and in what order. In practice, this will be used by MediumSystem
        // so that we can easily construct the \c mixNumberDensv argument of updateGasState.
        int h;
        string grainType;

        // Representative sizes. Will be used naively in the gas code (processes are calculated
        // separately for each size given here, and then summed). Integrating everything over the
        // grain size distribution is not doable for most processes in the gas code.
        Array sizev;
        Array numberDensRatiov;     // Number of grains 'per hydrogen atom' for each size
        std::vector<Array> qabsvv;  // Q_abs(a, nu), indexed on (size, frequency)
    };

    /** Initialize the gas module. Should be called exactly once, before any other functions of
        this class are used. The first argument should be the grid used for meanIntensityv, i.e.
        the radiation field wavelength grid. The second argument describes the dust population
        properties. The wavelength grid given to the third argument will be used to calculate the
        emission. */
    static void initialize(const Array& lambdav, const std::vector<DustInfo>&, const Array& emissionWLG);

    /** Function which should be called exactly one, before the code finishes. After calling
        finalize(), all info about the gas will be lost. */
    static void finalize();

    /** This function allocates space for the resuls of updateGasState to be stored. The gas is
        initialized as transparent. If updateGasState has not yet been called for a certain \c m,
        the returned values will be zero. */
    static void allocateGasStates(size_t num);

    /** This function returns \c true if the grain type described by the given string (get it from
        MultiGrainDustMic::populationGrainType(c)), is supported by the gas module. */
    static bool hasGrainTypeSupport(const string& populationGrainType);

    /** This function returns the \c h indices that were contained in the dust info objects at
        initialization, in the right order. They can be used to address the right density for the
        mixNumberDensv argument of updateGasState(). */
    static const vector<int>& hIndices();

    /** This function prepares the arguments for \c GasInterface::updateGasState() and calls it.
        The \c m argument indicates which of the allocated gas states should be updated. The number
        density of the gas medium and the radiation field (mean intensity) should be given as
        arguments, as well as the number densities for the dust media. The dust number densities
        for each size will be calculated by multiplying the numbers in \c mixNumberDensv with the
        \c numberDensRatiov of each \c DustInfo. This argument should contain the hydrogen number
        density factor \c n for medium \c h of cell \c m. This function is thread safe as long as
        each thread works on different \c m. */
    static void updateGasState(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv);

    /** This function re-reruns the calculation of the gas equilibrium for a cell, and returns many
        diagnostic values as a vector. The exact contents of this vector are subject to frequent
        change. Use \c diagnosticNames() to find out what these values mean. */
    static vector<double> diagnostics(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv);

    /** This functions returns a vector containing the names of the diagnostics returned by \c
        diagnostics(). */
    static vector<string> diagnosticNames();

    /** This function synchronizes the gas properties that need to be present at each process. It
        should be called by each process after they have finished working on updateGasState in
        parallel. TODO: communicate gas state contents. */
    static void communicateResults();

    /** This function initializes all values of the gas properties to zero. The function should be
        called before updateGasState() is called using multiprocessing. By starting from a clean
        slate each time, the communicatin in \c communicateResults() can simply be done using a
        sum. If some of the current results are needed to calculate the new state, they be copied
        over to a different variable; something like that should be implemented here. */
    static void clearResults();

    /** This function returns the gas temperature for cell m that resulted from \c updateGasState(m). */
    static double temperature(int m);

    /** This function returns the density of H+ in cell m that resulted from \c updateGasState(m). */
    static double np(int m);

    /** This function returns the density of H atoms in cell m that resulted from \c updateGasState(m). */
    static double nH(int m);

    /** This function returns the density of H2 molecules in cell m that resulted from \c updateGasState(m). */
    static double nH2(int m);

    /** This function calculates the opacity of the gas at a given wavelength for state \c m. If \c
        updateGasState() has not yet been called for cell \c m, the return value will be 0. Behind
        the scenes, the opacity is implemented using a table, of which the rows are filled when \c
        updateGasState is called. */
    static double opacityAbs(double lambda, int m);

    /** This function is the same as \c opacityAbs(double, int), but with a known wavelength index
        (which can be calculated using \c indexForLambda(). This provides a more efficient way of
        calculating the opacity for many cells at the same wavelength. */
    static double opacityAbs(int ell, int m);

    /** This function returns the index in the internal opacity table for the given wavelength, so
        that the overload of \c opacityAbs() below can be used. */
    static int indexForLambda(double lambda);

    /** This function calculates the emissivity (W m-3 sr-1 m-1) on the wavelength grid that was
        given at initialization, using the gas state stored at index m. This quantity is evaluated
        for the wavelengths in the emissivity wavelength grid given at initialization. The result
        only includes the continuum emission, and hence only the basic details of the gas are
        needed. */
    static Array emissivity(int m);

    /** This function calculates the emissivity (W m-3 sr-1 m-1) on the wavelength grid that was
        given at initialization, with emission lines included. Some quantities have to be
        recalculated to obtain the line ratios. Therefore the radiation field and grain densities
        have to be provided in the same format as updateGasState */
    static Array emissivity(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv);
};

#endif
