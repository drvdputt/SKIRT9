#include "Gas.hpp"
#include "Array.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"
#include "GrainInterface.hpp"
#include "NR.hpp"
#include "ProcessManager.hpp"
#include "StringUtils.hpp"
#include "Table.hpp"

#ifdef BUILD_WITH_GAS
#    include "GasDiagnostics.hpp"
#    include "GasInterface.hpp"
#    include <chrono>
#    include <iostream>
#    include <vector>
#endif

// note: (at least) all code using stuff related to GasInterface or the GasState vector should be
// wrapped in ifdef BUILD_WITH_GAS. The rest can be compiled and remain uninitialized. (Leaving the
// rest be might help a bit with the proliferation of macro's here.

namespace
{
    // set during initialize
    // ---------------------

    Array _lambdav;   // wavelengths given at initialization
    Array _olambdav;  // wavelengths for determining the index in the opacity table
    Array _elambdav;  // wavelengths for calculating the emission
#ifdef BUILD_WITH_GAS
    GasModule::GasInterface* _gi;  // instance of GasInterface
#endif
    std::vector<Gas::DustInfo> _dustinfov;  // information about the dust populations in the simulation
    vector<int> _hIndices;                  // indices which can be returned to help MediumSystem
    int _ip{-1}, _iH{-1}, _iH2{-1};         // indices to retrieve densities from gas states

// set per cell during updateGasState
// ----------------------------------
#ifdef BUILD_WITH_GAS
    std::vector<GasModule::GasState> _statev;  // result of the equilibrium calculation for each cell
#endif
    Table<2> _opacityvv;  // opacity(m, ell) for each cell m and wavelength ell

// utility functions
// -----------------

// translate a SKIRT grain type into one of the built-in grain types
#ifdef BUILD_WITH_GAS
    GasModule::GrainTypeLabel stringToGrainTypeLabel(const string& populationGrainType)
    {
        if (StringUtils::contains(populationGrainType, "Silicate"))
            return GasModule::GrainTypeLabel::SIL;
        else if (StringUtils::contains(populationGrainType, "Graphite")
                 || StringUtils::contains(populationGrainType, "PAH"))
            return GasModule::GrainTypeLabel::CAR;
        else
            return GasModule::GrainTypeLabel::OTHER;
    }
#endif

    // dlambda = - (c / nu^2) dnu
    // dnu = - (c / lambda^2) dlambda
    //
    // jnu dnu = jnu (- c / lambda^2) dlambda = jlambda dlambda
    // --> jnu = lambda^2 / c * jlambda
    // and conversion is involution (its own inverse)

    // convert from a quantity per x to a quantity per (c x^-1)
    Array x_to_cxm1(const Array& xv, const Array& quantity_xv)
    {
        int numx = xv.size();
        Array quantity_cxm1v(numx);
        for (int ix = 0; ix < numx; ix++)
        {
            double x = xv[ix];
            int icxm1 = numx - 1 - ix;
            quantity_cxm1v[icxm1] = x * x / Constants::c() * quantity_xv[ix];
        }
        return quantity_cxm1v;
    }

    Array lambdaToNu(const Array& lambdav, const Array& quantityPerLambda)
    {
        return x_to_cxm1(lambdav, quantityPerLambda);
    }

    Array nuToLambda(const Array& nuv, const Array& quantityPerNu) { return x_to_cxm1(nuv, quantityPerNu); }

    /** Convert array of x to array of 1 / x, with the elements ordered the other way round */
    Array invertAndFlip(const Array& xv)
    {
        int size = xv.size();
        Array xv_inv_flip(size);
        for (int i = 0; i < size; i++) xv_inv_flip[i] = 1. / xv[size - 1 - i];
        return xv_inv_flip;
    }

    // Get an array of grain number densities [cm-3] corresponding to dustInfo i with mix number
    // density mixNumberDens [m-3] (n in MediumSystem)
    Array mixNumberDensToGrainDensityv(int i, double mixNumberDens)
    {
        return _dustinfov[i].numberDensRatiov * mixNumberDens * 1.e-6;
    }

    // thread locals for efficiency
    // ----------------------------

    // properly initialized and modified by the first call to setThreadLocalGrainDensities
#ifdef BUILD_WITH_GAS
    thread_local GasModule::GrainInterface t_grainInterface;
    thread_local bool t_gr_is_ready{false};
    void setThreadLocalGrainDensities(const Array& mixNumberDensv, bool verbose)
    {
        // initialize when a thread meets this function for the first time (i.e. no populations are
        // present yet)
        if (!t_gr_is_ready)
        {
            for (size_t i = 0; i < _dustinfov.size(); i++)
            {
                // Just use 30 as the initial guess for the dust temperature, since SKIRT doesn't really
                // support calculating the dust temperature for individual sizes.
                Array temperaturev(30., _dustinfov[i].sizev.size());
                // Set the grain number densities using the number density of the mix (fictional H
                // density), and change unit from m-3 to cm-3
                Array densityv = mixNumberDensToGrainDensityv(i, mixNumberDensv[i]);
                t_grainInterface.addPopulation(stringToGrainTypeLabel(_dustinfov[i].grainType), _dustinfov[i].sizev,
                                               densityv, temperaturev, _gi->iFrequencyv(), _dustinfov[i].qabsvv);
                t_gr_is_ready = true;
            }
        }
        else
        {
            // simply change the number densities of the populations added in the block above
            for (size_t i = 0; i < _dustinfov.size(); i++)
            {
                const Array& densityv = mixNumberDensToGrainDensityv(i, mixNumberDensv[i]);
                if (verbose)
                {
                    std::cout << "pop " << i << " grain sizes:";
                    for (double d : _dustinfov[i].sizev) std::cout << ' ' << d;
                    std::cout << "\npop " << i << " grain densities:";
                    for (double d : densityv) std::cout << ' ' << d;
                    std::cout << '\n';
                }
                t_grainInterface.changePopulationDensityv(i, densityv);
            }
        }
    }
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::initialize(const Array& lambdav, const std::vector<DustInfo>& dustinfov, const Array& emissionWLG)
{
#ifdef BUILD_WITH_GAS
    if (_gi) FATALERROR("Gas module should be initialized exactly once");

    _lambdav = lambdav;
    _elambdav = emissionWLG;
    _dustinfov = dustinfov;
    _hIndices.reserve(_dustinfov.size());

    // Change the units of the dust properties from SI to cgs
    for (Gas::DustInfo& d : _dustinfov)
    {
        _hIndices.push_back(d.h);
        // Change size unit m to cm
        d.sizev *= 100.;
        // Flip the qabs arrays, because frequencies. This happens in-place.
        for (size_t b = 0; b < d.qabsvv.size(); b++) std::reverse(std::begin(d.qabsvv[b]), std::end(d.qabsvv[b]));
    }

    // Calculate the input radiation field / output opacity frequency grid
    Array iFrequencyv = Constants::c() * invertAndFlip(_lambdav);

    // Calculate the output emissivity frequency grid
    Array eFrequencyv = Constants::c() * invertAndFlip(_elambdav);

    // derive a wavelength grid that will be used for converting a wavelength to an index in the
    // opacity table (copied from DustMix)
    int numLambda = _lambdav.size();
    _olambdav.resize(numLambda);
    _olambdav[0] = lambdav[0];
    for (int ell = 1; ell != numLambda; ++ell) _olambdav[ell] = sqrt(lambdav[ell] * lambdav[ell - 1]);

    // Turn off error handling (otherwise, gas module can call abort)
    GasModule::GasInterface::errorHandlersOff();
    // Initialize the gas module
    _gi = new GasModule::GasInterface(iFrequencyv, iFrequencyv, eFrequencyv);

    // retrieve some useful indices
    _ip = _gi->index("H+");
    _iH = _gi->index("H");
    _iH2 = _gi->index("H2");
#else
    throw FATALERROR("SKIRT was built without gas support!")(void) lambdav;
    (void)dustinfov;
    (void)emissionWLG;
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::finalize()
{
#ifdef BUILD_WITH_GAS
    delete _gi;
    _gi = nullptr;
#endif
}

////////////////////////////////////////////////////////////////////

void Gas::allocateGasStates(size_t num)
{
#ifdef BUILD_WITH_GAS
    _statev.resize(num);
    _opacityvv.resize(num, _gi->oFrequencyv().size());
#else
    (void)num;
#endif
}

////////////////////////////////////////////////////////////////////

bool Gas::hasGrainTypeSupport(const string& populationGrainType)
{
#ifdef BUILD_WITH_GAS
    return stringToGrainTypeLabel(populationGrainType) != GasModule::GrainTypeLabel::OTHER;
#else
    (void)populationGrainType;
    return false;
#endif
}

////////////////////////////////////////////////////////////////////

const vector<int>& Gas::hIndices()
{
    return _hIndices;
}

////////////////////////////////////////////////////////////////////

#ifdef BUILD_WITH_GAS
namespace
{
    // implementation of updateGasState, with optional GasDiagnostics pointer
    void updateGasState_impl(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv,
                             GasModule::GasDiagnostics* gasDiagnostics)
    {

        auto start = std::chrono::high_resolution_clock::now();
        const Array& iFrequencyv = _gi->iFrequencyv();

        if (iFrequencyv.size() != meanIntensityv.size())
            throw FATALERROR("Something went wrong with the wavelength/frequency grids");

        Array jnu = lambdaToNu(_lambdav, meanIntensityv);
        // unit conversion:
        // for gas module: erg s-1 cm-2 sr-1 Hz-1
        // for skirt     : J   s-1 m-2  sr-1 Hz-1
        //                 7   0   -4
        jnu *= 1.e3;

        size_t countzeros = 0;
        for (size_t i = 0; i < iFrequencyv.size(); i++)
            if (jnu[i] <= 0) countzeros++;

        bool verbose = !(m % 300);
        if (verbose && countzeros) std::cout << countzeros << " zeros in cell " << m << '\n';

        // prepare grain info for this cell
        setThreadLocalGrainDensities(mixNumberDensv, verbose);

        // calculate the equilibrium
        _gi->updateGasState(_statev[m], n * 1.e-6, jnu, t_grainInterface, gasDiagnostics);

        // calculate and store the opacity; the opacity table is indexed on wavelength, so we need to
        // flip the result around
        const Array& opacity_nu = _gi->opacityWithLines(_statev[m], jnu, t_grainInterface, true, false, true);
        for (size_t ell = 0; ell < opacity_nu.size(); ell++)
            _opacityvv(m, ell) = opacity_nu[opacity_nu.size() - 1 - ell];

        if (verbose)
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            std::cout << "gas sample " << m << " n " << n * 1.e-6 << " time " << duration.count() << " ms.\n";
            std::cout << _gi->quickInfo(_statev[m], jnu) << '\n';
        }
    }
}
#endif

////////////////////////////////////////////////////////////////////

void Gas::updateGasState(int m, double n, const Array& meanIntensityv, const Array& mixNumberDensv)
{
#ifdef BUILD_WITH_GAS
    updateGasState_impl(m, n, meanIntensityv, mixNumberDensv, nullptr);
#else
    (void)m;
    (void)n;
    (void)meanIntensityv;
    (void)mixNumberDensv;
#endif
}

////////////////////////////////////////////////////////////////////

vector<double> Gas::diagnostics(int m, double n, const Array& meanintensityv, const Array& mixNumberDensv)
{
    vector<double> result;
#ifdef BUILD_WITH_GAS
    // recalculate the gas state and extract diagnostics (expensive)
    GasModule::GasDiagnostics gasDiagnostics;
    updateGasState_impl(m, n, meanintensityv, mixNumberDensv, &gasDiagnostics);

    // Gather the results; note that each map (e.g. gd.heating() and gd.cooling()) should contain
    // the contributions in the same order each time. I don't know if this is guaranteed by the c++
    // spec, but if it isn't, then I will need some other way to make sure that diagnosticNames()
    // and diagnostics() match.

    // heating, cooling, reaction rates
    result.reserve(20);
    for (auto& pair : gasDiagnostics.heating()) result.emplace_back(pair.second);
    for (auto& pair : gasDiagnostics.cooling()) result.emplace_back(pair.second);
    for (auto& d : gasDiagnostics.reactionRates()) result.emplace_back(d);
#else
    (void)m;
    (void)n;
    (void)meanIntensityv;
    (void)mixNumberDensv;
#endif
    return result;
}

////////////////////////////////////////////////////////////////////

vector<string> Gas::diagnosticNames()
{
    vector<string> result;
#ifdef BUILD_WITH_GAS
    // do a dummy calculation to get a gasdiagnostics object and figure out the names
    GasModule::GasState gasState;
    GasModule::GasDiagnostics gasDiagnostics;
    GasModule::GrainInterface grainInterface;
    _gi->updateGasState(gasState, 0, Array(_gi->iFrequencyv().size()), grainInterface, &gasDiagnostics);

    // heating, cooling, reaction rates, hopefully in the same order as
    result.reserve(20);
    for (auto& pair : gasDiagnostics.heating()) result.emplace_back(pair.first);
    for (auto& pair : gasDiagnostics.cooling()) result.emplace_back(pair.first);
    for (auto& s : gasDiagnostics.reactionNames()) result.emplace_back(s);
#endif
    return result;
}

////////////////////////////////////////////////////////////////////

void Gas::communicateResults()
{
    ProcessManager::sumToAll(_opacityvv.data());
}

////////////////////////////////////////////////////////////////////

void Gas::clearResults()
{
    _opacityvv.setToZero();
}

////////////////////////////////////////////////////////////////////

double Gas::temperature(int m)
{
#ifdef BUILD_WITH_GAS
    return _statev[m].temperature();
#else
    (void)m;
    return 0;
#endif
}

////////////////////////////////////////////////////////////////////

namespace
{
    // avoid some duplication of the ifdef stuff for the functions below
    double density(int m, int index)
    {
#ifdef BUILD_WITH_GAS
        return _statev[m].density(index);
#else
        (void)m;
        (void)index;
        return 0;
#endif
    }
}

////////////////////////////////////////////////////////////////////

double Gas::np(int m)
{
    return density(m, _ip);
}

////////////////////////////////////////////////////////////////////

double Gas::nH(int m)
{
    return density(m, _iH);
}

////////////////////////////////////////////////////////////////////

double Gas::nH2(int m)
{
    return density(m, _iH2);
}

////////////////////////////////////////////////////////////////////

double Gas::opacityAbs(double lambda, int m)
{
    return _opacityvv(m, indexForLambda(lambda));
}

////////////////////////////////////////////////////////////////////

double Gas::opacityAbs(int ell, int m)
{
    return _opacityvv(m, ell);
}

////////////////////////////////////////////////////////////////////

int Gas::indexForLambda(double lambda)
{
    return NR::locateClip(_olambdav, lambda);
}

////////////////////////////////////////////////////////////////////

Array Gas::emissivity(int m)
{
#ifdef BUILD_WITH_GAS
    return nuToLambda(_gi->eFrequencyv(), _gi->emissivityBasic(_statev[m], true));
#else
    (void)m;
    return Array;
#endif
}

////////////////////////////////////////////////////////////////////
