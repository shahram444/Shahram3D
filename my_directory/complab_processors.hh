#include "../defineKinetics.hh"
#include "../neuralNetwork.hh"
#include <random>
using namespace plb;
typedef double T;

#define DESCRIPTOR descriptors::D3Q19Descriptor  // Cs2 = 1/3
#define BGK descriptors::AdvectionDiffusionD3Q7Descriptor // Cs2 = 1/4

#define thrd DBL_EPSILON

/* ===============================================================================================================
   =============================================== DATA PROCESSORS ===============================================
   =============================================================================================================== */

   // nonlocal biomass reaction

   // status : the generic status of the current basic solution
   // 1 (GLP_UNDEF)  - solution is undefined
   // 2 (GLP_FEAS)   - solution is feasible
   // 3 (GLP_INFEAS) - solution is infeasible
   // 4 (GLP_NOFEAS) - problem has no feasible solution
   // 5 (GLP_OPT)    - solution is optimal
   // 6 (GLP_UNBND)  - problem has unbounded solution

   // fmin : the current values of the objective function : growth rate
   // time : time consumed for run_glpk()
   // mem  : memory used for run_glpk()

   // std::vector<double> xmin(ncol), lambda(nrow), redcosts(ncol);
   // xmin : the primal value vector of the structural variable (ncol) : flux vector
   // lambda : the dual value vector (i.e. reduced cost) of the auxiliary variable (nrow)
   // redcosts : the dual value vector (i.e. reduced cost) of the structural variable (ncol).

template<typename T, template<typename U> class Descriptor>
class run_FBA : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    // run_FBA( plint nx_, plint subsNum_, plint bioNum_, T dt_, T max_bMassRho_, std::vector<T> vec1_mu_, std::vector< std::vector<T> > vec2_maxUptake_, std::vector< std::vector<T> > vec2_maxRelease_, std::vector<bool> vec1_fixLB_, std::vector<bool> vec1_fixC_, std::vector< std::vector<T> > vec2_Kc_, std::vector< std::vector<int> > vec2_subsLoc_, std::vector<int> vec1_nrow_,
    //              std::vector<int> vec1_ncol_, std::vector<glp_prob *> vec1_lp_, std::vector<int> method_, std::vector<int> isMIP_, std::vector<glp_smcp> sParam_, std::vector<glp_iocp> iParam_, std::vector<plint> objLoc_, std::vector<std::vector<T>> vec2_Vmax_, std::vector<plint> pore_, plint solid_, plint bb_ )
    // : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), dt(dt_), max_bMassRho(max_bMassRho_), vec1_mu(vec1_mu_), vec2_maxUptake(vec2_maxUptake_), vec2_maxRelease(vec2_maxRelease_), vec1_fixLB(vec1_fixLB_), vec1_fixC(vec1_fixC_), vec2_Kc(vec2_Kc_), vec2_subsLoc(vec2_subsLoc_), vec1_nrow(vec1_nrow_),vec1_ncol(vec1_ncol_),
    //   vec1_lp(vec1_lp_), method(method_), isMIP(isMIP_), sParam(sParam_), iParam(iParam_), objLoc(objLoc_), vec2_Vmax(vec2_Vmax_), pore(pore_), solid(solid_), bb(bb_)
    run_FBA(plint nx_,  plint subsNum_, plint bioNum_, T dt_, T max_bMassRho_, std::vector<T> vec1_mu_, std::vector< std::vector<T> > vec2_maxUptake_, std::vector< std::vector<T> > vec2_maxRelease_, std::vector<bool> vec1_fixLB_, std::vector<bool> vec1_fixC_, std::vector< std::vector<T> > vec2_Kc_, std::vector< std::vector<int> > vec2_subsLoc_, std::vector<int> vec1_nrow_,
        std::vector<int> vec1_ncol_, std::vector<glp_prob*> vec1_lp_, std::vector<int> method_, std::vector<int> isMIP_, std::vector<glp_smcp> sParam_, std::vector<glp_iocp> iParam_, std::vector<plint> objLoc_, std::vector<std::vector<T>> vec2_Vmax_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), dt(dt_), max_bMassRho(max_bMassRho_), vec1_mu(vec1_mu_), vec2_maxUptake(vec2_maxUptake_), vec2_maxRelease(vec2_maxRelease_), vec1_fixLB(vec1_fixLB_), vec1_fixC(vec1_fixC_), vec2_Kc(vec2_Kc_), vec2_subsLoc(vec2_subsLoc_), vec1_nrow(vec1_nrow_), vec1_ncol(vec1_ncol_),
        vec1_lp(vec1_lp_), method(method_), isMIP(isMIP_), sParam(sParam_), iParam(iParam_), objLoc(objLoc_), vec2_Vmax(vec2_Vmax_)
    {}
    // substrate lattices and then glpk_bio lattices. bio-lattices contain domain material numbers
    // dt in seconds
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        // T thrd = 1e-12; // T thrd = DBL_EPSILON;
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            if (absX > 0 && absX < nx - 1) {
                std::vector<Dot3D> vec_offset;
                plint maskLloc = subsNum + bioNum;
                for (plint iT = 0; iT < maskLloc + 1; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT])); }
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {

                        // for (plint iS=0; iS<subsNum; ++iS) {
                        //     if (lattices[iS]->get(iX+vec_offset[iS].x,iY+vec_offset[iS].y).computeDensity() < 0) {
                        //         Array<T,5> g;
                        //         g[0]=(T) (-1)/3; g[1]=g[2]=g[3]=g[4]=(T) (-1)/6;
                        //         lattices[iS]->get(iX+vec_offset[iS].x,iY+vec_offset[iS].y).setPopulations(g);
                        //     }
                        // }

                        // check if the current location is the grid cell with bmass
                        T bioT = 0.;
                        plint num_of_microbes = 0;
                        std::vector<T> bmass;
                        std::vector<plint> bmLattice;
                        for (plint iM = 0; iM < bioNum; ++iM) {
                            T bval = lattices[subsNum + iM]->get(iX + vec_offset[subsNum + iM].x, iY + vec_offset[subsNum + iM].y, iZ + vec_offset[subsNum + iM].z).computeDensity();
                            if (bval > thrd) {
                                ++num_of_microbes; bioT += bval; // unit: kgDW/m3
                                bmass.push_back(bval); // vector storing biomass (unit: kgDW/m3)
                                bmLattice.push_back(iM); // vector storing the lattice order of bmass vector
                            }
                        }
                        if (num_of_microbes > 0) {
                            std::vector<T> conc;
                            // construct the concentration vector
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                T c0 = lattices[iS]->get(iX + vec_offset[iS].x, iY + vec_offset[iS].y, iZ + vec_offset[iS].z).computeDensity(); // mmol/L
                                if (c0 <= thrd) { c0 = 0; }
                                conc.push_back(c0); // mmol/L
                            }

                            bool errflag = 1, skipflag = 0;
                            plint loopcount = 0;
                            std::vector<bool> errchk(subsNum, 0);
                            std::vector<T> vec_dc(subsNum, 0.), vec1_grate(num_of_microbes);
                            std::vector< std::vector<bool> > tmpchk(num_of_microbes);
                            for (plint iter = 0; iter < num_of_microbes; ++iter) { tmpchk[iter] = std::vector<bool>(subsNum); }
                            while (errflag == 1) {
                                // execute run_glpk for each metabolic model (loop over each microbe)
                                std::vector<std::vector<T>> vec2_flux(num_of_microbes), vec2_lb;
                                for (plint iB = 0; iB < num_of_microbes; ++iB) {
                                    T g_rate = 0; plint bloc = bmLattice[iB];
                                    std::vector<T> vec1_lb, vec1_ub, vec1_flux;
                                    for (plint iS = 0; iS < subsNum; ++iS) {
                                        vec1_ub.push_back(vec2_maxRelease[bloc][iS]);
                                        if (vec2_subsLoc[bloc][iS] > 0) {
                                            T c2f = -conc[iS] / (conc[iS] + vec2_Kc[bloc][iS]);
                                            if (errchk[iS] == 0) {
                                                if (vec2_Vmax[bloc][iS] > thrd) { c2f *= vec2_Vmax[bloc][iS]; }
                                                else { c2f *= conc[iS] / bmass[iB] / dt * 3600; }
                                            }
                                            else {
                                                T common_biomass = 0;
                                                for (plint iB1 = 0; iB1 < num_of_microbes; ++iB1) {
                                                    if (tmpchk[iB1][iS] == 1) { common_biomass += bmass[iB1]; }
                                                }
                                                c2f *= conc[iS] / common_biomass / dt * 3600;
                                            }
                                            if (c2f < vec2_maxUptake[bloc][iS]) { c2f = vec2_maxUptake[bloc][iS]; }
                                            if (c2f > vec1_ub[iS]) { skipflag = 1; break; }
                                            if (vec1_fixLB[iS] == 1 && c2f < conc[iS]) { c2f = -conc[iS]; }
                                            vec1_lb.push_back(c2f);
                                        }
                                        else { vec1_lb.push_back(-0.); } // to match the size of the vec1_lb with vec2_subsLoc[bloc]
                                    }
                                    if (skipflag == 1) { break; }
                                    int status, lpsolver = 1;
                                    T time, mem;
                                    std::vector<T> fluxes(vec1_ncol[bloc]), lambda(vec1_nrow[bloc]), redcosts(vec1_ncol[bloc]);

                                    // if (lexID == 0) {
                                    int glpkerr = run_glpk(vec1_lp[bloc], method[bloc], isMIP[bloc], vec1_nrow[bloc], vec1_ncol[bloc], lpsolver, sParam[bloc], iParam[bloc],
                                        fluxes, g_rate, status, vec2_subsLoc[bloc], vec1_lb, vec1_ub, lambda, redcosts, time, mem);
                                    if (!glpkerr) {
                                        if (status == 180 || status == 5) {
                                            // solution found
                                            for (plint iS = 0; iS < subsNum; ++iS) {
                                                if (vec2_maxUptake[bmLattice[iB]][iS] > 0) {
                                                    vec1_flux.push_back(fluxes[vec2_subsLoc[bloc][iS]]);
                                                }
                                                else {
                                                    vec1_flux.push_back(0.);
                                                }
                                            }
                                        }
                                        else {
                                            g_rate = 0.;
                                            for (plint iS = 0; iS < subsNum; ++iS) {
                                                vec1_flux.push_back(0.);
                                            }
                                        }
                                    }
                                    else {
                                        // pcout << "ERROR in run_glpk. errnum = ";
                                        if (glpkerr - 100 == GLP_EBADB) {
                                            std::cout << "invalid base.\n";
                                            std::cout << "Unable to start the search, because the initial basis specified in the problem object is invalid:" <<
                                                " the number of basic (auxiliary and structural) variables is not the same as the number of rows in the problem object.\n";
                                        }
                                        else if (glpkerr - 100 == GLP_ESING) { std::cout << "singular matrix.\n"; }
                                        else if (glpkerr - 100 == GLP_ECOND) { std::cout << "ill-conditioned matrix.\n"; }
                                        else if (glpkerr - 100 == GLP_EBOUND) { std::cout << "incorrect bounds.\n"; }
                                        else if (glpkerr - 100 == GLP_EFAIL) { std::cout << "solver failure.\n"; }
                                        else if (glpkerr - 100 == GLP_EOBJLL) { std::cout << "lower limit reached.\n"; }
                                        else if (glpkerr - 100 == GLP_EOBJUL) { std::cout << "upper limit reached.\n"; }
                                        else if (glpkerr - 100 == GLP_EITLIM) { std::cout << "iterations limit exceeded.\n"; }
                                        else if (glpkerr - 100 == GLP_ETMLIM) { std::cout << "time limit exceeded.\n"; }
                                        else if (glpkerr - 100 == GLP_ENODFS) { std::cout << "no dual feasible solution.\n"; }
                                        else if (glpkerr - 100 == GLP_ENOPFS) { std::cout << "no primal feasible solution.\n"; }
                                        g_rate = 0.;
                                        for (plint iS = 0; iS < subsNum; ++iS) {
                                            vec1_flux.push_back(0.);
                                        }
                                    }
                                    vec2_flux[iB] = std::vector<T>(subsNum);
                                    for (plint iS = 0; iS < subsNum; ++iS) {
                                        vec2_flux[iB][iS] = vec1_flux[iS];
                                    }
                                    vec1_grate[iB] = g_rate;
                                    vec2_lb.push_back(vec1_lb);
                                }

                                if (skipflag == 0) {
                                    // error handling
                                    for (plint iS = 0; iS < subsNum; ++iS) {
                                        if (vec1_fixLB[iS] == 0 && vec1_fixC[iS] == 0) {
                                            T dc = 0;
                                            for (plint iB = 0; iB < num_of_microbes; ++iB) {
                                                // if (vec2_maxUptake[bloc][iS] > 0) {
                                                    // f2c = f2c*vec2_Kc[bloc][iS]/(vec2_Vmax[bloc][iS]-f2c);
                                                dc += vec2_flux[iB][iS] / 3600 * bmass[iB] * dt;
                                                // }
                                            }
                                            if (dc + conc[iS] < -thrd && conc[iS]>thrd) {
                                                std::cout << "Concentration of the substrate " << iS << " goes negative (" << dc + conc[iS] << "). Adjust lower bounds and rerun GLPK (" << loopcount << ").\n";
                                                for (plint iB = 0; iB < num_of_microbes; ++iB) {
                                                    for (plint iSS = 0; iSS < subsNum; ++iSS) {
                                                        if (vec2_flux[iB][iSS] < 0) { tmpchk[iB][iSS] = 1; }
                                                        else { tmpchk[iB][iSS] = 0; }
                                                    }
                                                }
                                                errchk[iS] = 1;
                                                if (loopcount > 1) exit(EXIT_FAILURE);
                                            }
                                            else { vec_dc[iS] = dc; errchk[iS] = 0; }
                                        }
                                        else { vec_dc[iS] = 0; errchk[iS] = 0; }
                                    }
                                    plint errcount = 0;
                                    for (plint iS = 0; iS < subsNum; ++iS) {
                                        if (errchk[iS] == 0) { ++errcount; }
                                    }
                                    if (errcount == subsNum) { errflag = 0; }
                                    ++loopcount;
                                }
                                else { break; }
                            }

                            if (skipflag == 0) {
                                // update concentration
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    if (vec1_fixC[iS] == 0 && vec1_fixLB[iS] == 0) {
                                        T dc = vec_dc[iS];
                                        if (dc > thrd || dc < -thrd) {
                                            Array<T, 7> g;
                                            lattices[iS]->get(iX + vec_offset[iS].x, iY + vec_offset[iS].y, iZ + vec_offset[iS].z).getPopulations(g);
                                            g[0] += (T)(dc) / 4; g[1] += (T)(dc) / 8; g[2] += (T)(dc) / 8; g[3] += (T)(dc) / 8; g[4] += (T)(dc) / 8; g[5] += (T)(dc) / 8; g[6] += (T)(dc) / 8;
                                            lattices[iS]->get(iX + vec_offset[iS].x, iY + vec_offset[iS].y, iZ + vec_offset[iS].z).setPopulations(g);
                                        }
                                    }
                                }

                                // update biomass
                                for (plint iB = 0; iB < num_of_microbes; ++iB) {
                                    plint bloc = bmLattice[iB];
                                    T g_rate = vec1_grate[iB]; // unit: gDW/gDW/hr
                                    if (g_rate > -thrd && g_rate < thrd) {
                                        g_rate = -vec1_mu[bloc];
                                    }
                                    T dB = g_rate * dt / 3600 * bmass[iB]; // unit: kgDW/m3 (=gDW/L)
                                    if (bmass[iB] + dB <= 0) { dB = -bmass[iB]; }
                                    if (bmass[iB] + dB > max_bMassRho) { dB = max_bMassRho - bmass[iB]; }
                                    if (dB > thrd || dB < -thrd) {
                                        Array<T, 7> g;
                                        lattices[subsNum + bloc]->get(iX + vec_offset[subsNum + bloc].x, iY + vec_offset[subsNum + bloc].y, iZ + vec_offset[subsNum + bloc].z).getPopulations(g);
                                        g[0] += (T)(dB) / 4; g[1] += (T)(dB) / 8; g[2] += (T)(dB) / 8; g[3] += (T)(dB) / 8; g[4] += (T)(dB) / 8; g[5] += (T)(dB) / 8; g[6] += (T)(dB) / 8;
                                        lattices[subsNum + bloc]->get(iX + vec_offset[subsNum + bloc].x, iY + vec_offset[subsNum + bloc].y, iZ + vec_offset[subsNum + bloc].z).setPopulations(g);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual run_FBA<T, Descriptor>* clone() const {
        return new run_FBA<T, Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iT = 0; iT < subsNum + bioNum; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }
private:
    plint nx, subsNum, bioNum;
    T dt;
    T max_bMassRho;
    std::vector<T> vec1_mu;
    std::vector<std::vector<T>> vec2_maxUptake, vec2_maxRelease;
    std::vector<bool> vec1_fixLB;
    std::vector<bool> vec1_fixC;
    std::vector<std::vector<T>> vec2_Kc;
    std::vector<std::vector<int>> vec2_subsLoc;
    std::vector<int> vec1_nrow, vec1_ncol;
    std::vector<glp_prob*> vec1_lp;
    std::vector<int> method, isMIP;
    std::vector<glp_smcp> sParam;
    std::vector<glp_iocp> iParam;
    std::vector<plint> objLoc;
    std::vector<std::vector<T>> vec2_Vmax;
    // std::vector<plint> pore;
    // plint solid, bb;
};

template<typename T, template<typename U> class Descriptor>
class run_kinetics : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    run_kinetics(plint nx_, plint subsNum_, plint bioNum_, T dt_, std::vector<std::vector<T>> vec2_Kc_kns_, std::vector<T> vec1_mu_kns_, plint solid_, plint bb_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), dt(dt_), vec2_Kc_kns(vec2_Kc_kns_), vec1_mu_kns(vec1_mu_kns_), solid(solid_), bb(bb_)
    {}
    // substrate lattices and then bio-lattices. the mask number lattice is always the last.
    // dt in seconds, dx in meters
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            if (absX > 0 && absX < nx - 1) {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                        plint maskLloc = subsNum + bioNum;
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(iX + maskOffset.x, iY + maskOffset.y, iZ + maskOffset.z).computeDensity());
                        if (mask != solid && mask != bb) {
                            std::vector<Dot3D> vec_offset;
                            for (plint iT = 0; iT < maskLloc; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT])); }
                            std::vector<T> bmass, conc, subs_rate(subsNum, 0), bio_rate(bioNum, 0);
                            // construct the concentration vector
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                plint iXs = iX + vec_offset[iS].x, iYs = iY + vec_offset[iS].y, iZs = iZ + vec_offset[iS].z;
                                T c0 = lattices[iS]->get(iXs, iYs, iZs).computeDensity(); // mmol/L
                                if (c0 < thrd) { c0 = 0; }
                                conc.push_back(c0); // mmol/L
                            }
                            // construct the biomass vector
                            for (plint iM = 0; iM < bioNum; ++iM) {
                                plint iXb = iX + vec_offset[subsNum + iM].x, iYb = iY + vec_offset[subsNum + iM].y, iZb = iZ + vec_offset[subsNum + iM].z;
                                T b0 = lattices[subsNum + iM]->get(iXb, iYb, iZb).computeDensity(); // kgDW/m3 = gDW/L
                                if (b0 < thrd) { b0 = 0; }
                                bmass.push_back(b0); // vector storing biomass (unit: kgDW/m3)
                            }

                            defineRxnKinetics(bmass, conc, subs_rate, bio_rate, mask);

                            // update concentration
                            for (plint iS = 0; iS < subsNum; ++iS) {
                                // forward-Euler method
                                T dC = subs_rate[iS] * dt;
                                if (dC > thrd || dC < -thrd) {
                                    if (conc[iS] + dC <= 0) { dC = -conc[iS]; }
                                    Array<T, 7> g;
                                    plint iXt = iX + vec_offset[iS].x, iYt = iY + vec_offset[iS].y, iZt = iZ + vec_offset[iS].z;
                                    lattices[iS]->get(iXt, iYt, iZt).getPopulations(g);
                                    g[0] += (T)(dC) / 4; g[1] += (T)(dC) / 8; g[2] += (T)(dC) / 8; g[3] += (T)(dC) / 8; g[4] += (T)(dC) / 8; g[5] += (T)(dC) / 8; g[6] += (T)(dC) / 8;
                                    lattices[iS]->get(iXt, iYt, iZt).setPopulations(g);
                                }
                            }
                            // update biomass
                            for (plint iB = 0; iB < bioNum; ++iB) {
                                // forward-Euler method
                                T dB = bio_rate[iB] * dt; // unit: kgDW/m3 (=gDW/L)
                                if (dB > thrd || dB < -thrd) {
                                    if (bmass[iB] + dB <= 0) { dB = -bmass[iB]; }
                                    Array<T, 7> g;
                                    plint iXt = iX + vec_offset[subsNum + iB].x, iYt = iY + vec_offset[subsNum + iB].y, iZt = iZ + vec_offset[subsNum + iB].z;
                                    lattices[subsNum + iB]->get(iXt, iYt, iZt).getPopulations(g);
                                    g[0] += (T)(dB) / 4; g[1] += (T)(dB) / 8; g[2] += (T)(dB) / 8; g[3] += (T)(dB) / 8; g[4] += (T)(dB) / 8; g[5] += (T)(dB) / 8; g[6] += (T)(dB) / 8;
                                    lattices[subsNum + iB]->get(iXt, iYt, iZt).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual run_kinetics<T, Descriptor>* clone() const {
        return new run_kinetics<T, Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iT = 0; iT < subsNum + bioNum; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }
private:
    plint nx, subsNum, bioNum;
    T dt;
    std::vector<std::vector<T>> vec2_Kc_kns;
    std::vector<T> vec1_mu_kns;
    plint solid, bb;
};

template<typename T, template<typename U> class Descriptor>
class run_NN : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    run_NN(plint nx_, plint subsNum_, plint bioNum_, plint solid_, plint bb_, T dt_, std::vector<plint> reaction_type_, std::vector<std::vector<T>> vec2_Vmax_, std::vector<std::vector<T>> vec2_Kc_, std::vector<T> vec1_mu_, bool flag_)
        : nx(nx_), subsNum(subsNum_), bioNum(bioNum_), solid(solid_), bb(bb_), dt(dt_), reaction_type(reaction_type_), vec2_Vmax(vec2_Vmax_), vec2_Kc(vec2_Kc_), vec1_mu(vec1_mu_), flag(flag_)
    {}
    // substrate lattices and then bio-lattices. the mask number lattice is always the last.
    // dt in seconds, dx in meters
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            if (absX > 0 && absX < nx - 1) {
                for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                    for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                        plint maskLloc = subsNum + bioNum;
                        Dot3D maskOffset = computeRelativeDisplacement(*lattices[0], *lattices[maskLloc]);
                        plint mask = util::roundToInt(lattices[maskLloc]->get(iX + maskOffset.x, iY + maskOffset.y, iZ + maskOffset.z).computeDensity());
                        if (mask != solid && mask != bb) {
                            std::vector<plint> bmLoc; std::vector<T> bmass; std::vector<Dot3D> vec_offset;
                            for (plint iT = 0; iT < maskLloc; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT])); }
                            // construct the biomass vector
                            for (plint iM = 0; iM < bioNum; ++iM) {
                                plint iXb = iX + vec_offset[subsNum + iM].x, iYb = iY + vec_offset[subsNum + iM].y, iZb = iZ + vec_offset[subsNum + iM].z;
                                T b0 = lattices[subsNum + iM]->get(iXb, iYb, iZb).computeDensity(); // kgDW/m3 = gDW/L
                                if (b0 > thrd) { bmass.push_back(b0); bmLoc.push_back(iM); } // vector storing biomass (unit: kgDW/m3)
                            }
                            if (bmass.size() > 0) {
                                // construct the concentration vector
                                std::vector<T> conc;
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    plint iXs = iX + vec_offset[iS].x, iYs = iY + vec_offset[iS].y, iZs = iZ + vec_offset[iS].z;
                                    T c0 = lattices[iS]->get(iXs, iYs, iZs).computeDensity(); // mmol/L
                                    if (c0 < thrd) { c0 = 0; }
                                    conc.push_back(c0); // vector storing biomass (unit: kgDW/m3)
                                }

                                // construct the flux vectors
                                std::vector<std::vector<T>> Fin, Fout;
                                for (size_t iB = 0; iB < bmLoc.size(); ++iB) {
                                    plint iM = bmLoc[iB];
                                    std::vector<T> tmp;
                                    T b0 = bmass[iB];
                                    for (plint iS = 0; iS < subsNum; ++iS) {
                                        T c2f = conc[iS] / (conc[iS] + vec2_Kc[iM][iS]);
                                        if (vec2_Vmax[iM][iS] > thrd) {
                                            if (flag == 1) {
                                                T dc = c2f * vec2_Vmax[iM][iS] / 3600 * bmass[iB] * dt;
                                                if (conc[iS] - dc < thrd) { c2f *= conc[iS] / b0 * (3600 / dt); }
                                                else { c2f *= vec2_Vmax[iM][iS]; }
                                            }
                                            else { c2f *= vec2_Vmax[iM][iS]; }
                                        }
                                        else {
                                            if (b0 < thrd) { c2f = 0; }
                                            else { c2f *= conc[iS] / b0 * (3600 / dt); }
                                        }
                                        tmp.push_back(c2f); // vector storing substrate influx (unit: kgDW/m3)
                                    }
                                    Fin.push_back(tmp);
                                }
                                std::vector<T> bio_rate(bmass.size(), 0);
                                Fout = Fin;

                                defineNeuralNetwork(Fin, Fout, bio_rate);

                                std::vector<T> vec_dc;
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    // flux to conc
                                    T dc = 0;
                                    for (size_t iB = 0; iB < Fout.size(); ++iB) { dc += Fout[iB][iS] / 3600 * bmass[iB] * dt; }
                                    dc *= -1;
                                    if (dc + conc[iS] < -thrd) {
                                        std::cout << "Concentration of the substrate " << iS << " goes negative (" << dc + conc[iS] << "). Terminating the simulation.\n";
                                        exit(EXIT_FAILURE);
                                    }
                                    // update concentration
                                    if (dc > thrd || dc < -thrd) {
                                        Array<T, 7> g;
                                        lattices[iS]->get(iX + vec_offset[iS].x, iY + vec_offset[iS].y, iZ + vec_offset[iS].z).getPopulations(g);
                                        g[0] += (T)(dc) / 4; g[1] += (T)(dc) / 8; g[2] += (T)(dc) / 8; g[3] += (T)(dc) / 8; g[4] += (T)(dc) / 8; g[5] += (T)(dc) / 8; g[6] += (T)(dc) / 8;
                                        lattices[iS]->get(iX + vec_offset[iS].x, iY + vec_offset[iS].y, iZ + vec_offset[iS].z).setPopulations(g);
                                    }
                                }

                                // update biomass
                                for (size_t iB = 0; iB < bmLoc.size(); ++iB) {
                                    plint iM = bmLoc[iB];
                                    T g_rate = bio_rate[iB]; // unit: gDW/gDW/hr
                                    if (g_rate > -thrd && g_rate < thrd) { g_rate = -vec1_mu[iM]; }
                                    T dB = g_rate * dt / 3600 * bmass[iB]; // unit: kgDW/m3 (=gDW/L)
                                    if (bmass[iB] + dB <= 0) { dB = -bmass[iB]; }
                                    if (dB > thrd || dB < -thrd) {
                                        Array<T, 7> g;
                                        lattices[subsNum + iM]->get(iX + vec_offset[subsNum + iM].x, iY + vec_offset[subsNum + iM].y, iZ + vec_offset[subsNum + iM].z).getPopulations(g);
                                        g[0] += (T)(dB) / 4; g[1] += (T)(dB) / 8; g[2] += (T)(dB) / 8; g[3] += (T)(dB) / 8; g[4] += (T)(dB) / 8; g[5] += (T)(dB) / 8; g[6] += (T)(dB) / 8;
                                        lattices[subsNum + iM]->get(iX + vec_offset[subsNum + iM].x, iY + vec_offset[subsNum + iM].y, iZ + vec_offset[subsNum + iM].z).setPopulations(g);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual run_NN<T, Descriptor>* clone() const {
        return new run_NN<T, Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iT = 0; iT < subsNum + bioNum; ++iT) {
            modified[iT] = modif::staticVariables;
        }
    }
private:
    plint nx, subsNum, bioNum, solid, bb;
    T dt;
    std::vector<plint> reaction_type;
    std::vector<std::vector<T>> vec2_Vmax, vec2_Kc;
    std::vector<T> vec1_mu;
    bool flag;
};


// redistribute the excessive biomass (push)
template<typename T, template<typename U> class Descriptor>
class pushExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    pushExcessBiomass3D(T Bmax_, plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_, plint solid_, plint bb_, std::vector<plint> pore_)
        : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_), solid(solid_), bb(bb_), pore(pore_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-3] = total biomass lattice
    // lattices[len-2] = mask lattice
    // lattices[len-1] = age lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        plint distLloc = length - 1, maskLloc = length - 2, bMtLloc = length - 3, numbM = (length - 3) / 2;
        for (plint iL = 0; iL < length; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL])); }
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
            plint iXm = iX0 + vec_offset[maskLloc].x;
            for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                plint iYm = iY0 + vec_offset[maskLloc].y;
                for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                    plint iZm = iZ0 + vec_offset[maskLloc].z;
                    plint mask = util::roundToInt(lattices[maskLloc]->get(iXm, iYm, iZm).computeDensity());
                    if (mask != bb && mask != solid) {
                        plint iXt = iX0 + vec_offset[bMtLloc].x, iYt = iY0 + vec_offset[bMtLloc].y, iZt = iZ0 + vec_offset[bMtLloc].z;
                        T bMt = lattices[bMtLloc]->get(iXt, iYt, iZt).computeDensity();
                        if (bMt > Bmax) {
                            T bMd = bMt - Bmax;
                            if (bMd > thrd) {
                                plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                                std::vector<plint> delXYZ; plint nbrs = 0;
                                if (absX == bdryGap && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {
                                    nbrs = 3; // nbrs is the number of neighbors depending on the current location
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);
                                   

                                }
                                else if (absX == (nx - 1 - bdryGap) && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX > bdryGap && absX < (nx - 1 - bdryGap) && absY == 0 && absZ == 0) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX > bdryGap && absX < (nx - 1 - bdryGap) && absY == (ny - 1) && absZ == (nz - 1)) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == bdryGap && absY == 0 && absZ == 0) {
                                    nbrs = 2;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX == bdryGap && absY == (ny - 1) && absZ == (nz - 1)) {
                                    nbrs = 2;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == (nx - 1 - bdryGap) && absY == 0 && absZ == 0) {
                                    nbrs = 2;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX == (nx - 1 - bdryGap) && absY == ny - 1 && absZ == nz - 1) {
                                    nbrs = 2;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else {
                                    nbrs = 4;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                
                                // check if iX or iY is located on the subdomain boundaries (MPI)
                                plint bdry_dir = 0;
                                // if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY > 0 && absY < (ny-1) ) {
                                if (absX >= bdryGap && absX <= (nx - 1 - bdryGap)) {
                                    if (iX0 == domain.x1 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 1; }
                                    else if (iY0 == domain.y1 && iX0 > domain.x0 && iX0 < domain.x1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 2; }

                                    else if (iZ0 == domain.z1 && iX0 > domain.x0 && iX0 < domain.x1 && iY0 > domain.x0 && iY0 < domain.y1) { bdry_dir = 3; }

                                    else if (iX0 == domain.x0 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 4; }
                                    else if (iY0 == domain.y0 && iX0 > domain.x0 && iX0 < domain.x1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 5; }

                                    else if (iZ0 == domain.z0 && iX0 > domain.x0 && iX0 < domain.x1 && iY0 > domain.y0 && iY0 < domain.y1) { bdry_dir = 6; }

                                    else if (iX0 == domain.x1 && iY0 == domain.y1 && iZ0 == domain.z1) { bdry_dir = 7; }
                                    else if (iX0 == domain.x0 && iY0 == domain.y1 && iZ0 == domain.z1) { bdry_dir = 8; }
                                    else if (iX0 == domain.x0 && iY0 == domain.y0 && iZ0 == domain.z0) { bdry_dir = 9; }
                                    else if (iX0 == domain.x1 && iY0 == domain.y0 && iZ0 == domain.z0) { bdry_dir = 10; }



                                }








                                // search neighbors for redistribution
                                std::vector<plint> nbrsLocMask; plint nbrTlen = 0;
                                for (plint iT = 0; iT < nbrs; ++iT) {
                                    plint delx = delXYZ[iT * 3], dely = delXYZ[iT * 3 + 1], delz = delXYZ[iT * 3 + 1];
                                    plint nbrmask = util::roundToInt(lattices[maskLloc]->get(iXm + delx, iYm + dely, iZm + delz).computeDensity());
                                    // exclude wall boundaries
                                    if (nbrmask != bb && nbrmask != solid) {
                                        ++nbrTlen;
                                        nbrsLocMask.push_back(delx); nbrsLocMask.push_back(dely); nbrsLocMask.push_back(dely); nbrsLocMask.push_back(nbrmask);
                                    }
                                }
                                // shuffle a number vector to randomly select a neighbor
                                std::vector<plint> randArray;
                                if (nbrTlen > 1) {
                                    for (plint iR = 0; iR < nbrTlen; ++iR) { randArray.push_back(iR); }
                                    std::random_shuffle(randArray.begin(), randArray.end());
                                }
                                // first trial to redistribute excess biomass
                                bool chk = 0;
                                if (nbrTlen > 0) {
                                    for (plint iT = 0; iT < nbrTlen; ++iT) {
                                        // select a neighboring grid cell
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = randArray[iT]; }
                                        plint delx = nbrsLocMask[3 * randLoc], dely = nbrsLocMask[3 * randLoc + 1], delz = nbrsLocMask[3 * randLoc + 1], nbrmask = nbrsLocMask[3 * randLoc + 2];
                                        T nbrbMt = 0;
                                        bool nbrmaskflag = 0;
                                        if (nbrmask == bb || nbrmask == solid) { nbrmaskflag = 1; }
                                        else { for (size_t iP = 0; iP < pore.size(); ++iP) { if (nbrmask == pore[iP]) { nbrmaskflag = 1; break; } } }
                                        if (nbrmaskflag == 0) { nbrbMt = lattices[bMtLloc]->get(iXt + delx, iYt + dely, iZt + delz).computeDensity(); }
                                        if (nbrbMt < Bmax) { // else, select a new neighbor
                                            // redefine the push direction indicator based on the direction of the chosen neighbor
                                            plint push_dir = 0;
                                            if (bdry_dir > 0) {
                                                if (delx == 1 && (bdry_dir == 1 || bdry_dir == 7 || bdry_dir == 10)) { push_dir = 1; }
                                                else if (dely == 1 && (bdry_dir == 2 || bdry_dir == 7 || bdry_dir == 8)) { push_dir = 2; }
                                                else if (delz == 1 && (bdry_dir == 3 || bdry_dir == 7 || bdry_dir == 8)) { push_dir = 3; }
                                                else if (delx == -1 && (bdry_dir == 4 || bdry_dir == 8 || bdry_dir == 9)) { push_dir = 4; }
                                                else if (dely == -1 && (bdry_dir == 5 || bdry_dir == 9 || bdry_dir == 10)) { push_dir = 5; }
                                                else if (delz == -1 && (bdry_dir == 6 || bdry_dir == 9 || bdry_dir == 10)) { push_dir = 6; }

                                                else { push_dir = 0; }
                                            }
                                            T hold_capacity = Bmax - nbrbMt;
                                            T partial_bMassT = bMd;
                                            if (bMd > hold_capacity) { partial_bMassT = hold_capacity; bMd -= hold_capacity; }
                                            else { bMd = 0; chk = 1; }
                                            for (plint iM = 0; iM < numbM; ++iM) {
                                                Array<T, 7> g;
                                                plint iXb1 = iX0 + vec_offset[iM].x, iYb1 = iY0 + vec_offset[iM].y, iZb1 = iZ0 + vec_offset[iM].z; // original biomass lattice
                                                T shove_bmass = (lattices[iM]->get(iXb1, iYb1, iZb1).computeDensity()) / bMt * partial_bMassT;
                                                if (shove_bmass > thrd) {
                                                    lattices[iM]->get(iXb1, iYb1, iZb1).getPopulations(g);
                                                    g[0] -= (T)(shove_bmass) / 4; g[1] -= (T)(shove_bmass) / 8; g[2] -= (T)(shove_bmass) / 8; g[3] -= (T)(shove_bmass) / 8;  g[4] -= (T)(shove_bmass) / 8; g[5] -= (T)(shove_bmass) / 8;   g[6] -= (T)(shove_bmass) / 8;
                                                    lattices[iM]->get(iXb1, iYb1, iZb1).setPopulations(g);
                                                    // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                                    if (push_dir > 0) {
                                                        plint iXb2 = iX0 + vec_offset[iM + numbM].x; plint iYb2 = iY0 + vec_offset[iM + numbM].y; plint iZb2 = iZ0 + vec_offset[iM + numbM].z;
                                                        g[0] = (T)(shove_bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(shove_bmass - 1) / 8;
                                                        lattices[iM + numbM]->get(iXb2, iYb2, iZb2).setPopulations(g);
                                                        lattices[iM + numbM]->get(iXb2, iYb2, iZb2).getDynamics().setOmega(push_dir);
                                                    }
                                                    // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                                    else {
                                                        lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).getPopulations(g);
                                                        g[0] += (T)(shove_bmass) / 4; g[1] += (T)(shove_bmass) / 8; g[2] += (T)(shove_bmass) / 8; g[3] += (T)(shove_bmass) / 8; g[4] += (T)(shove_bmass) / 8; g[5] += (T)(shove_bmass) / 8; g[6] += (T)(shove_bmass) / 8;
                                                        lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).setPopulations(g);
                                                    }
                                                }
                                            }
                                            bMt -= partial_bMassT;
                                        }
                                        // break out of the for loop if the excess biomass has been redistributed
                                        if (chk == 1) { break; }
                                    }
                                }
                                else {
                                    std::cout << "there is no neighbor for biomass redistribution. terminating the simulation.\n";
                                    exit(EXIT_FAILURE);
                                }
                                // redistribute the remaining biomass (this is the most time consuming part)
                                if (chk == 0) {
                                    /* purely random version (too slow)*/
                                    // plint randLoc=0;
                                    // if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                    // plint delx = nbrsLocMask[randLoc*3], dely = nbrsLocMask[randLoc*3+1];

                                    /* use the age lattice */

                                    std::vector<T> tmp1Vector;
                                    plint tmp1Len = 0;
                                    plint iXd = iX0 + vec_offset[distLloc].x, iYd = iY0 + vec_offset[distLloc].y, iZd = iZ0 + vec_offset[distLloc].z;
                                    plint id0 = util::roundToInt(lattices[distLloc]->get(iXd, iYd, iZd).computeDensity());
                                    for (plint iT = 0; iT < nbrTlen; ++iT) {
                                        plint delx = nbrsLocMask[iT * 3], dely = nbrsLocMask[iT * 3 + 1], delz = nbrsLocMask[iT * 3 + 1];
                                        plint id1 = util::roundToInt(lattices[distLloc]->get(iXd + delx, iYd + dely, iZd + delz).computeDensity());
                                        if (id0 > id1) { tmp1Vector.push_back(delx); tmp1Vector.push_back(dely); tmp1Vector.push_back(delz); ++tmp1Len; }
                                    }
                                    plint delx, dely, delz, push_dir = 0;
                                    if (tmp1Len > 1) {
                                        plint randLoc = std::rand() % tmp1Len;
                                        delx = tmp1Vector[randLoc * 2]; dely = tmp1Vector[randLoc * 2 + 1]; delz = tmp1Vector[randLoc * 2 + 1];
                                    }
                                    else if (tmp1Len == 1) {
                                        delx = tmp1Vector[0]; dely = tmp1Vector[1]; delz = tmp1Vector[1];
                                    }
                                    else { // tmp1Len == 0, purely random redistribution
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                        delx = nbrsLocMask[randLoc * 3]; dely = nbrsLocMask[randLoc * 3 + 1]; delz = nbrsLocMask[randLoc * 3 + 1];
                                    }
                                    if (bdry_dir > 0) {
                                        if (delx == 1 && (bdry_dir == 1 || bdry_dir == 7 || bdry_dir == 10)) { push_dir = 1; }
                                        else if (dely == 1 && (bdry_dir == 2 || bdry_dir == 7 || bdry_dir == 8)) { push_dir = 2; }
                                        else if (delz == 1 && (bdry_dir == 3 || bdry_dir == 7 || bdry_dir == 8)) { push_dir = 3; }
                                        else if (delx == -1 && (bdry_dir == 4 || bdry_dir == 8 || bdry_dir == 9)) { push_dir = 4; }
                                        else if (dely == -1 && (bdry_dir == 5 || bdry_dir == 9 || bdry_dir == 10)) { push_dir = 5; }
                                        else if (delz == -1 && (bdry_dir == 6 || bdry_dir == 9 || bdry_dir == 10)) { push_dir = 6; }

                                        else { push_dir = 0; }
                                    }


                                    for (plint iM = 0; iM < numbM; ++iM) {
                                        Array<T, 7> g;
                                        plint iXb1 = iX0 + vec_offset[iM].x, iYb1 = iY0 + vec_offset[iM].y, iZb1 = iZ0 + vec_offset[iM].z;; // original biomass lattice
                                        T shove_bmass = (lattices[iM]->get(iXb1, iYb1, iZb1).computeDensity()) / bMt * bMd;
                                        lattices[iM]->get(iXb1, iYb1, iZb1).getPopulations(g);
                                        g[0] -= (T)(shove_bmass) / 4; g[1] -= (T)(shove_bmass) / 8; g[2] -= (T)(shove_bmass) / 8; g[3] -= (T)(shove_bmass) /8 ; g[4] -= (T)(shove_bmass) / 8; g[5] -= (T)(shove_bmass) / 8; g[6] -= (T)(shove_bmass) / 8;
                                        lattices[iM]->get(iXb1, iYb1, iZb1).setPopulations(g);
                                        // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                        if (push_dir > 0) {
                                            plint iXb2 = iX0 + vec_offset[iM + numbM].x; plint iYb2 = iY0 + vec_offset[iM + numbM].y; plint iZb2 = iZ0 + vec_offset[iM + numbM].z;
                                            g[0] = (T)(shove_bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(shove_bmass - 1) / 8;
                                            lattices[iM + numbM]->get(iXb2, iYb2, iZb2).setPopulations(g);
                                            lattices[iM + numbM]->get(iXb2, iYb2, iZb2).getDynamics().setOmega(push_dir);
                                        }
                                        // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                        else {
                                            lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).getPopulations(g);
                                            g[0] += (T)(shove_bmass) / 4; g[1] += (T)(shove_bmass) / 8; g[2] += (T)(shove_bmass) / 8; g[3] += (T)(shove_bmass) / 8; g[4] += (T)(shove_bmass) / 8; g[5] += (T)(shove_bmass) / 8; g[6] += (T)(shove_bmass) / 8;
                                            lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).setPopulations(g);
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
    }
        virtual BlockDomain::DomainT appliesTo() const {
            // Don't apply to envelope, because nearest neighbors need to be accessed.
            return BlockDomain::bulk;
        }
        virtual pushExcessBiomass3D<T, Descriptor>* clone() const {
            return new pushExcessBiomass3D<T, Descriptor>(*this);
        }
        void getTypeOfModification(std::vector<modif::ModifT>&modified) const {
            plint numbM = (length - 3) / 2;
            for (plint iB = 0; iB < numbM; ++iB) {
                modified[iB] = modif::staticVariables;
                modified[iB + numbM] = modif::allVariables;
            }
            modified[length - 3] = modif::nothing;
            modified[length - 2] = modif::nothing;
            modified[length - 1] = modif::nothing;
        }
private:
       T Bmax;
       plint nx, ny, nz, bdryGap, length, solid, bb;
       std::vector<plint> pore;
};



//new version
// redistribute the excessive biomass (push)
template<typename T, template<typename U> class Descriptor>
class pushExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    pushExcessBiomass3D(T Bmax_, plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_, plint solid_, plint bb_, std::vector<plint> pore_)
    : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_), solid(solid_), bb(bb_), pore(pore_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-4)] = copy biomass lattices
    // lattices[len-4] = total biomass lattice
    // lattices[len-3] = mask lattice
    // lattices[len-2] = age lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        plint distLloc = length-1, maskLloc = length-2, bMtLloc = length-3, numbM = (length-3)/2;
        for (plint iL=0; iL<length; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0],*lattices[iL])); }
        Dot3D absoluteOffset = lattices[0]->getLocation();
        for (plint iX0=domain.x0; iX0<=domain.x1; ++iX0) {
            plint iXm = iX0+vec_offset[maskLloc].x;
            for (plint iY0=domain.y0; iY0<=domain.y1; ++iY0) {
                plint iYm = iY0+vec_offset[maskLloc].y;
                for (plint iZ0=domain.z0; iZ0<=domain.z1; ++iZ0) {
                    plint iZm = iZ0+vec_offset[maskLloc].z;
                    plint mask = util::roundToInt( lattices[maskLloc]->get(iXm,iYm,iZm).computeDensity() );
                    if (mask != bb && mask != solid) {
                        plint iXt = iX0 + vec_offset[bMtLloc].x, iYt = iY0 + vec_offset[bMtLloc].y, iZt = iZ0 + vec_offset[bMtLloc].z;
                        T bMt = lattices[bMtLloc]->get(iXt,iYt,iZt).computeDensity();
                        if (bMt > Bmax) {
                            T bMd = bMt-Bmax;
                            if ( bMd > thrd ) {
                                plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                                std::vector<plint> delXYZ; plint nbrs = 0;
                                if ( absX == bdryGap && absY > 0 && absY < (ny-1) && absZ > 0 && absZ < (nz-1) ) {
                                    nbrs = 5; // nbrs is the number of neighbors depending on the current location
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    else if ( absX == (nx-1-bdryGap) && absY > 0 && absY < (ny-1) && absZ > 0 && absZ < (nz-1) ) {
                                    nbrs = 5;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    else if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY == 0 && absZ > 0 && absZ < (nz-1) ) {
                                    nbrs = 5;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    else if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY == (ny-1) && absZ > 0 && absZ < (nz-1) ) {
                                    nbrs = 5;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    else if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY > 0 && absY < (ny-1) && absZ == 0 ) {
                                    nbrs = 5;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    else if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY > 0 && absY < (ny-1) && absZ == (nz-1) ) {
                                    nbrs = 5;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    else if ( absX == bdryGap && absY == 0 && absZ == 0 ) {
                                    nbrs = 3;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == bdryGap && absY == (ny-1) && absZ == 0 ) {
                                    nbrs = 3;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == bdryGap && absY == 0 && absZ == (nz-1) ) {
                                    nbrs = 3;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == bdryGap && absY == (ny-1) && absZ == (nz-1) ) {
                                    nbrs = 3;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == (nx-1-bdryGap) && absY == 0 && absZ == 0 ) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == (nx-1-bdryGap) && absY == (ny-1) && absZ == 0 ) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == (nx-1-bdryGap) && absY == 0 && absZ == (nz-1) ) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0);
                                    }
                                    else if ( absX == (nx-1-bdryGap) && absY == (ny-1) && absZ == (nz-1) ) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    }
                                    else {
                                    nbrs = 6;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0);
                                    delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                    }
                                    plint bdry_dir = 0;
                                    // if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY > 0 && absY < (ny-1) ) {
                                    if ( absX >= bdryGap && absX <= (nx-1-bdryGap) ) {
                                        if (iX0 == domain.x1 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 1; }
                                        else if (iY0 == domain.y1 && iX0 > domain.x0 && iX0 < domain.x1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 2; }
                                        else if (iX0 == domain.x0 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 3; }
                                        else if (iY0 == domain.y0 && iX0 > domain.x0 && iX0 < domain.x1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 4; }
                                        else if (iX0 == domain.x1 && iY0 == domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 5; }
                                        else if (iX0 == domain.x0 && iY0 == domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 6; }
                                        else if (iX0 == domain.x0 && iY0 == domain.y0 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 7; }
                                        else if (iX0 == domain.x1 && iY0 == domain.y0 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 8; }
                                        else if (iX0 > domain.x0 && iX0 < domain.x1 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 == domain.z1) { bdry_dir = 9; }
                                        else if (iX0 > domain.x0 && iX0 < domain.x1 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 == domain.z0) { bdry_dir = 10; }
                                    }
                                    std::vector<plint> nbrsLocMask; plint nbrTlen = 0;
                                    // search neighbors for redistribution
                                    std::vector<plint> nbrsLocMask; plint nbrTlen = 0;
                                    for (plint iT=0; iT<nbrs; ++iT) {
                                        plint delx = delXYZ[iT*3], dely = delXYZ[iT*3+1], delz = delXYZ[iT*3+2];
                                        plint nbrmask = util::roundToInt( lattices[maskLloc]->get(iXm+delx,iYm+dely,iZm+delz).computeDensity() );
                                        // exclude wall boundaries
                                        if (nbrmask != bb && nbrmask != solid) {
                                            ++nbrTlen;
                                            nbrsLocMask.push_back(delx); nbrsLocMask.push_back(dely); nbrsLocMask.push_back(delz); nbrsLocMask.push_back(nbrmask);
                                        }
                                    }
                                    // shuffle a number vector to randomly select a neighbor
                                    std::vector<plint> randArray;
                                    if (nbrTlen > 1) {
                                        for (plint iR=0; iR<nbrTlen; ++iR) { randArray.push_back(iR); }
                                        std::random_shuffle(randArray.begin(), randArray.end());
                                    }
                                     // first trial to redistribute excess biomass
                                bool chk = 0;
                                if (nbrTlen > 0) {
                                    for (plint iT = 0; iT < nbrTlen; ++iT) {
                                        // select a neighboring grid cell
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = randArray[iT]; }
                                        plint delx = nbrsLocMask[3 * randLoc], dely = nbrsLocMask[3 * randLoc + 1], delz = nbrsLocMask[3 * randLoc + 1], nbrmask = nbrsLocMask[3 * randLoc + 2];
                                        T nbrbMt = 0;
                                        bool nbrmaskflag = 0;
                                        if (nbrmask == bb || nbrmask == solid) { nbrmaskflag = 1; }
                                        else { for (size_t iP = 0; iP < pore.size(); ++iP) { if (nbrmask == pore[iP]) { nbrmaskflag = 1; break; } } }
                                        if (nbrmaskflag == 0) { nbrbMt = lattices[bMtLloc]->get(iXt + delx, iYt + dely, iZt + delz).computeDensity(); }
                                        if (nbrbMt < Bmax) { // else, select a new neighbor
                                            // redefine the push direction indicator based on the direction of the chosen neighbor
                                           plint push_dir = 0;
                                        if (bdry_dir > 0) {
                                            if ( delx == 1 && (bdry_dir == 1 || bdry_dir == 5 || bdry_dir == 8) && delz == 0) { push_dir = 1; }
                                            else if ( dely == 1 && (bdry_dir == 2 || bdry_dir == 5 || bdry_dir == 6) && delz == 0) { push_dir = 2; }
                                            else if ( delx == -1 && (bdry_dir == 3 || bdry_dir == 6 || bdry_dir == 7) && delz == 0) { push_dir = 3; }
                                            else if ( dely == -1 && (bdry_dir == 4 || bdry_dir == 7 || bdry_dir == 8) && delz == 0) { push_dir = 4; }
                                            else if ( delz == 1 && (bdry_dir == 9 || bdry_dir == 10 || bdry_dir == 11 || bdry_dir == 12)) { push_dir = 5; }
                                            else if ( delz == -1 && (bdry_dir == 13 || bdry_dir == 14 || bdry_dir == 15 || bdry_dir == 16)) { push_dir = 6; }
                                            else { push_dir = 0; }
                                        }

                                            T hold_capacity = Bmax - nbrbMt;
                                            T partial_bMassT = bMd;
                                            if (bMd > hold_capacity) { partial_bMassT = hold_capacity; bMd -= hold_capacity; }
                                            else { bMd = 0; chk = 1; }
                                            for (plint iM = 0; iM < numbM; ++iM) {
                                                Array<T, 7> g;
                                                plint iXb1 = iX0 + vec_offset[iM].x, iYb1 = iY0 + vec_offset[iM].y, iZb1 = iZ0 + vec_offset[iM].z; // original biomass lattice
                                                T shove_bmass = (lattices[iM]->get(iXb1, iYb1, iZb1).computeDensity()) / bMt * partial_bMassT;
                                                if (shove_bmass > thrd) {
                                                    lattices[iM]->get(iXb1, iYb1, iZb1).getPopulations(g);
                                                    g[0] -= (T)(shove_bmass) / 4; g[1] -= (T)(shove_bmass) / 8; g[2] -= (T)(shove_bmass) / 8; g[3] -= (T)(shove_bmass) / 8;  g[4] -= (T)(shove_bmass) / 8; g[5] -= (T)(shove_bmass) / 8;   g[6] -= (T)(shove_bmass) / 8;
                                                    lattices[iM]->get(iXb1, iYb1, iZb1).setPopulations(g);
                                                    // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                                    if (push_dir > 0) {
                                                        plint iXb2 = iX0 + vec_offset[iM + numbM].x; plint iYb2 = iY0 + vec_offset[iM + numbM].y; plint iZb2 = iZ0 + vec_offset[iM + numbM].z;
                                                        g[0] = (T)(shove_bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(shove_bmass - 1) / 8;
                                                        lattices[iM + numbM]->get(iXb2, iYb2, iZb2).setPopulations(g);
                                                        lattices[iM + numbM]->get(iXb2, iYb2, iZb2).getDynamics().setOmega(push_dir);
                                                    }
                                                    // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                                    else {
                                                        lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).getPopulations(g);
                                                        g[0] += (T)(shove_bmass) / 4; g[1] += (T)(shove_bmass) / 8; g[2] += (T)(shove_bmass) / 8; g[3] += (T)(shove_bmass) / 8; g[4] += (T)(shove_bmass) / 8; g[5] += (T)(shove_bmass) / 8; g[6] += (T)(shove_bmass) / 8;
                                                        lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).setPopulations(g);
                                                    }
                                                }
                                            }
                                            bMt -= partial_bMassT;
                                        }
                                        // break out of the for loop if the excess biomass has been redistributed
                                        if (chk == 1) { break; }
                                    }
                                }
                                else {
                                    std::cout << "there is no neighbor for biomass redistribution. terminating the simulation.\n";
                                    exit(EXIT_FAILURE);
                                }
                                // redistribute the remaining biomass (this is the most time consuming part)
                                if (chk == 0) {
                                    /* purely random version (too slow)*/
                                    // plint randLoc=0;
                                    // if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                    // plint delx = nbrsLocMask[randLoc*3], dely = nbrsLocMask[randLoc*3+1];

                                    /* use the age lattice */

                                    std::vector<T> tmp1Vector;
                                    plint tmp1Len = 0;
                                    plint iXd = iX0 + vec_offset[distLloc].x, iYd = iY0 + vec_offset[distLloc].y, iZd = iZ0 + vec_offset[distLloc].z;
                                    plint id0 = util::roundToInt(lattices[distLloc]->get(iXd, iYd, iZd).computeDensity());
                                    for (plint iT = 0; iT < nbrTlen; ++iT) {
                                        plint delx = nbrsLocMask[iT * 3], dely = nbrsLocMask[iT * 3 + 1], delz = nbrsLocMask[iT * 3 + 1];
                                        plint id1 = util::roundToInt(lattices[distLloc]->get(iXd + delx, iYd + dely, iZd + delz).computeDensity());
                                        if (id0 > id1) { tmp1Vector.push_back(delx); tmp1Vector.push_back(dely); tmp1Vector.push_back(delz); ++tmp1Len; }
                                    }
                                    plint delx, dely, delz, push_dir = 0;
                                    if (tmp1Len > 1) {
                                        plint randLoc = std::rand() % tmp1Len;
                                        delx = tmp1Vector[randLoc * 2]; dely = tmp1Vector[randLoc * 2 + 1]; delz = tmp1Vector[randLoc * 2 + 1];
                                    }
                                    else if (tmp1Len == 1) {
                                        delx = tmp1Vector[0]; dely = tmp1Vector[1]; delz = tmp1Vector[1];
                                    }
                                    else { // tmp1Len == 0, purely random redistribution
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                        delx = nbrsLocMask[randLoc * 3]; dely = nbrsLocMask[randLoc * 3 + 1]; delz = nbrsLocMask[randLoc * 3 + 1];
                                    }
                                    if (bdry_dir > 0) {
                                        if ( delx == 1 && (bdry_dir == 1 || bdry_dir == 5 || bdry_dir == 8) && delz == 0) { push_dir = 1; }
                                        else if ( dely == 1 && (bdry_dir == 2 || bdry_dir == 5 || bdry_dir == 6) && delz == 0) { push_dir = 2; }
                                        else if ( delx == -1 && (bdry_dir == 3 || bdry_dir == 6 || bdry_dir == 7) && delz == 0) { push_dir = 3; }
                                        else if ( dely == -1 && (bdry_dir == 4 || bdry_dir == 7 || bdry_dir == 8) && delz == 0) { push_dir = 4; }
                                        else if ( delz == 1 && (bdry_dir == 9 || bdry_dir == 13 || bdry_dir == 15 || bdry_dir == 17 || bdry_dir == 18 || bdry_dir == 20) ) { push_dir = 5; }
                                        else if ( delz == -1 && (bdry_dir == 10 || bdry_dir == 11 || bdry_dir == 14 || bdry_dir == 16 || bdry_dir == 19 || bdry_dir == 21) ) { push_dir = 6; }
                                        else { push_dir = 0; }
                                    }

                                    for (plint iM = 0; iM < numbM; ++iM) {
                                        Array<T, 7> g;
                                        plint iXb1 = iX0 + vec_offset[iM].x, iYb1 = iY0 + vec_offset[iM].y, iZb1 = iZ0 + vec_offset[iM].z;; // original biomass lattice
                                        T shove_bmass = (lattices[iM]->get(iXb1, iYb1, iZb1).computeDensity()) / bMt * bMd;
                                        lattices[iM]->get(iXb1, iYb1, iZb1).getPopulations(g);
                                        g[0] -= (T)(shove_bmass) / 4; g[1] -= (T)(shove_bmass) / 8; g[2] -= (T)(shove_bmass) / 8; g[3] -= (T)(shove_bmass) /8 ; g[4] -= (T)(shove_bmass) / 8; g[5] -= (T)(shove_bmass) / 8; g[6] -= (T)(shove_bmass) / 8;
                                        lattices[iM]->get(iXb1, iYb1, iZb1).setPopulations(g);
                                        // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                        if (push_dir > 0) {
                                            plint iXb2 = iX0 + vec_offset[iM + numbM].x; plint iYb2 = iY0 + vec_offset[iM + numbM].y; plint iZb2 = iZ0 + vec_offset[iM + numbM].z;
                                            g[0] = (T)(shove_bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(shove_bmass - 1) / 8;
                                            lattices[iM + numbM]->get(iXb2, iYb2, iZb2).setPopulations(g);
                                            lattices[iM + numbM]->get(iXb2, iYb2, iZb2).getDynamics().setOmega(push_dir);
                                        }
                                        // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                        else {
                                            lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).getPopulations(g);
                                            g[0] += (T)(shove_bmass) / 4; g[1] += (T)(shove_bmass) / 8; g[2] += (T)(shove_bmass) / 8; g[3] += (T)(shove_bmass) / 8; g[4] += (T)(shove_bmass) / 8; g[5] += (T)(shove_bmass) / 8; g[6] += (T)(shove_bmass) / 8;
                                            lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).setPopulations(g);
                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }
        }
    }
        virtual BlockDomain::DomainT appliesTo() const {
            // Don't apply to envelope, because nearest neighbors need to be accessed.
            return BlockDomain::bulk;
        }
        virtual pushExcessBiomass3D<T, Descriptor>* clone() const {
            return new pushExcessBiomass3D<T, Descriptor>(*this);
        }
        void getTypeOfModification(std::vector<modif::ModifT>&modified) const {
            plint numbM = (length - 3) / 2;
            for (plint iB = 0; iB < numbM; ++iB) {
                modified[iB] = modif::staticVariables;
                modified[iB + numbM] = modif::allVariables;
            }
            modified[length - 3] = modif::nothing;
            modified[length - 2] = modif::nothing;
            modified[length - 1] = modif::nothing;
        }
private:
       T Bmax;
       plint nx, ny, nz, bdryGap, length, solid, bb;
       std::vector<plint> pore;
};

    // redistribute the excessive biomass (push)
    template<typename T, template<typename U> class Descriptor>
    class halfPushExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
    {
    public:
        halfPushExcessBiomass3D(T Bmax_, plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_, plint solid_, plint bb_, std::vector<plint> pore_)
            : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_), solid(solid_), bb(bb_), pore(pore_)
        {}
        // lattices[0~(#ofbM-1)] = original biomass lattices
        // lattices[#ofbM~(len-3)] = copy biomass lattices
        // lattices[len-3] = total biomass lattice
        // lattices[len-2] = mask lattice
        // lattices[len-1] = dist lattice
        virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
            std::vector<Dot3D> vec_offset;
            plint distLloc = length - 1, maskLloc = length - 2, bMtLloc = length - 3, numbM = (length - 3) / 2;
            for (plint iL = 0; iL < length; ++iL) {
                vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL]));
            }
            Dot3D absoluteOffset = lattices[0]->getLocation();
            for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                plint iXm = iX0 + vec_offset[maskLloc].x;
                for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                    plint iYm = iY0 + vec_offset[maskLloc].y;
                    for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                        plint iZm = iZ0 + vec_offset[maskLloc].z;
                        plint mask = util::roundToInt(lattices[maskLloc]->get(iXm, iYm, iZm).computeDensity());
                        bool maskflag = 0;
                        if (mask == bb || mask == solid) { maskflag = 1; }
                        else {
                            for (size_t iP = 0; iP < pore.size(); ++iP) {
                                if (mask == pore[iP]) {
                                    maskflag = 1;
                                    break;
                                }
                            }
                        }
                        if (maskflag == 0) {
                            plint iXt = iX0 + vec_offset[bMtLloc].x, iYt = iY0 + vec_offset[bMtLloc].y, iZt = iZ0 + vec_offset[bMtLloc].z;
                            T bMt = lattices[bMtLloc]->get(iXt, iYt, iZt).computeDensity();
                            if (bMt > Bmax) {
                                T bMd = bMt * 0.5;
                                plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                                std::vector<plint> delXYZ; plint nbrs = 0;
                                if (absX == bdryGap && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {
                                    nbrs = 3; // nbrs is the number of neighbors depending on the current location
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);
                              

                                }
                                else if (absX == (nx - 1 - bdryGap) && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX > bdryGap && absX < (nx - 1 - bdryGap) && absY == 0 && absZ == 0) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX > bdryGap && absX < (nx - 1 - bdryGap) && absY == (ny - 1) && absZ == (nz - 1)) {
                                    nbrs = 3;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == bdryGap && absY == 0 && absZ == 0) {
                                    nbrs = 2;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX == bdryGap && absY == (ny - 1) && absZ == (nz - 1)) {
                                    nbrs = 2;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == (nx - 1 - bdryGap) && absY == 0 && absZ == 0) {
                                    nbrs = 2;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX == (nx - 1 - bdryGap) && absY == ny - 1 && absZ == nz - 1) {
                                    nbrs = 2;
                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else {
                                    nbrs = 4;
                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }

                                // check if iX or iY is located on the subdomain boundaries (MPI)
                                plint bdry_dir = 0;
                                // if ( absX > bdryGap && absX < (nx-1-bdryGap) && absY > 0 && absY < (ny-1) ) {
                                if (absX >= bdryGap && absX <= (nx - 1 - bdryGap)) {
                                    if (iX0 == domain.x1 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 1; }
                                    else if (iY0 == domain.y1 && iX0 > domain.x0 && iX0 < domain.x1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 2; }

                                    else if (iZ0 == domain.z1 && iX0 > domain.x0 && iX0 < domain.x1 && iY0 > domain.x0 && iY0 < domain.y1) { bdry_dir = 3; }

                                    else if (iX0 == domain.x0 && iY0 > domain.y0 && iY0 < domain.y1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 4; }
                                    else if (iY0 == domain.y0 && iX0 > domain.x0 && iX0 < domain.x1 && iZ0 > domain.z0 && iZ0 < domain.z1) { bdry_dir = 5; }

                                    else if (iZ0 == domain.z0 && iX0 > domain.x0 && iX0 < domain.x1 && iY0 > domain.y0 && iY0 < domain.y1) { bdry_dir = 6; }

                                    else if (iX0 == domain.x1 && iY0 == domain.y1 && iZ0 == domain.z1) { bdry_dir = 7; }
                                    else if (iX0 == domain.x0 && iY0 == domain.y1 && iZ0 == domain.z1) { bdry_dir = 8; }
                                    else if (iX0 == domain.x0 && iY0 == domain.y0 && iZ0 == domain.z0) { bdry_dir = 9; }
                                    else if (iX0 == domain.x1 && iY0 == domain.y0 && iZ0 == domain.z0) { bdry_dir = 10; }



                                }
                                // search neighbors for redistribution
                                std::vector<plint> nbrsLocMask; plint nbrTlen = 0;
                                for (plint iT = 0; iT < nbrs; ++iT) {
                                    plint delx = delXYZ[iT * 3], dely = delXYZ[iT * 3 + 1], delz = delXYZ[iT * 3 + 2];
                                    plint nbrmask = util::roundToInt(lattices[maskLloc]->get(iXm + delx, iYm + dely, iZm + delz).computeDensity());
                                    // exclude wall boundaries
                                    if (nbrmask != bb && nbrmask != solid) {
                                        ++nbrTlen;

                                        nbrsLocMask.push_back(delx); nbrsLocMask.push_back(dely); nbrsLocMask.push_back(delz); nbrsLocMask.push_back(nbrmask);
                                    }
                                }
                                // shuffle a number vector to randomly select a neighbor
                                std::vector<plint> randArray;
                                if (nbrTlen > 1) {
                                    for (plint iR = 0; iR < nbrTlen; ++iR) { randArray.push_back(iR); }
                                    std::random_shuffle(randArray.begin(), randArray.end());
                                }
                                // first trial to redistribute excess biomass
                                bool chk = 0;
                                if (nbrTlen > 0) {
                                    for (plint iT = 0; iT < nbrTlen; ++iT) {
                                        // select a neighboring grid cell
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = randArray[iT]; }
                                        plint delx = nbrsLocMask[3 * randLoc], dely = nbrsLocMask[3 * randLoc + 1], delz = nbrsLocMask[3 * randLoc + 2], nbrmask = nbrsLocMask[3 * randLoc + 3];
                                        T nbrbMt = 0; bool nbrmaskflag = 0;
                                        if (nbrmask == bb || nbrmask == solid) { nbrmaskflag = 1; }
                                        else { for (size_t iP = 0; iP < pore.size(); ++iP) { if (nbrmask == pore[iP]) { nbrmaskflag = 1; break; } } }
                                        if (nbrmaskflag == 0) { nbrbMt = lattices[bMtLloc]->get(iXt + delx, iYt + dely, iZt + delz).computeDensity(); }
                                        if (nbrbMt + bMd <= Bmax) { // else, select a new neighbor
                                            // redefine the push direction indicator based on the direction of the chosen neighbor
                                            plint push_dir = 0;
                                            if (bdry_dir > 0) {
                                                if (delx == 1 && (bdry_dir == 1 || bdry_dir == 7 || bdry_dir == 10)) { push_dir = 1; }
                                                else if (dely == 1 && (bdry_dir == 2 || bdry_dir == 7 || bdry_dir == 8)) { push_dir = 2; }
                                                else if (delz == 1 && (bdry_dir == 3 || bdry_dir == 7 || bdry_dir == 8)) { push_dir = 3; }
                                                else if (delx == -1 && (bdry_dir == 4 || bdry_dir == 8 || bdry_dir == 9)) { push_dir = 4; }
                                                else if (dely == -1 && (bdry_dir == 5 || bdry_dir == 9 || bdry_dir == 10)) { push_dir = 5; }
                                                else if (delz == -1 && (bdry_dir == 6 || bdry_dir == 9 || bdry_dir == 10)) { push_dir = 6; }

                                                else { push_dir = 0; }
                                            }
                                            for (plint iM = 0; iM < numbM; ++iM) {
                                                Array<T, 7> g;
                                                plint iXb1 = iX0 + vec_offset[iM].x, iYb1 = iY0 + vec_offset[iM].y, iZb1 = iZ0 + vec_offset[iM].z; // original biomass lattice
                                                T shove_bmass = (lattices[iM]->get(iXb1, iYb1, iZb1).computeDensity()) * 0.5;
                                                lattices[iM]->get(iXb1, iYb1, iZb1).getPopulations(g);
                                                g[0] -= (T)(shove_bmass) / 4; g[1] -= (T)(shove_bmass) / 8; g[2] -= (T)(shove_bmass) / 8; g[3] -= (T)(shove_bmass) / 6; g[4] -= (T)(shove_bmass) / 6; g[5] -= (T)(shove_bmass) / 8; g[6] -= (T)(shove_bmass) / 8;
                                                lattices[iM]->get(iXb1, iYb1, iZb1).setPopulations(g);
                                                // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                                if (push_dir > 0) {
                                                    plint iXb2 = iX0 + vec_offset[iM + numbM].x; plint iYb2 = iY0 + vec_offset[iM + numbM].y; plint iZb2 = iZ0 + vec_offset[iM + numbM].z;
                                                    g[0] = (T)(shove_bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(shove_bmass - 1) / 8;
                                                    lattices[iM + numbM]->get(iXb2, iYb2, iZb2).setPopulations(g);
                                                    lattices[iM + numbM]->get(iXb2, iYb2, iZb2).getDynamics().setOmega(push_dir);
                                                }
                                                // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                                else {
                                                    lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).getPopulations(g);
                                                    g[0] += (T)(shove_bmass) / 4; g[1] += (T)(shove_bmass) / 8; g[2] += (T)(shove_bmass) / 8; g[3] += (T)(shove_bmass) / 8; g[4] += (T)(shove_bmass) / 8; g[5] += (T)(shove_bmass) / 8; g[6] += (T)(shove_bmass) / 8;
                                                    lattices[iM]->get(iXb1 + delx, iYb1 + dely, iYb1 + delz).setPopulations(g);
                                                }
                                            }
                                        }
                                        // break out of the for loop if the excess biomass has been redistributed
                                        if (chk == 1) { break; }
                                    }
                                }
                                else {
                                    std::cout << "there is no neighbor for biomass redistribution. terminating the simulation.\n";
                                    exit(EXIT_FAILURE);
                                }
                                // redistribute the remaining biomass (this is the most time consuming part)
                                plint push_dir = 0, delx = 0, dely = 0, delz = 0;
                                if (chk == 0) {
                                    /* use the distance lattice */
                                    std::vector<T> tmp1Vector;
                                    plint tmp1Len = 0;
                                    plint iXd = iX0 + vec_offset[distLloc].x, iYd = iY0 + vec_offset[distLloc].y, iZd = iZ0 + vec_offset[distLloc].z;
                                    plint id0 = util::roundToInt(lattices[distLloc]->get(iXd, iYd, iZd).computeDensity());
                                    for (plint iT = 0; iT < nbrTlen; ++iT) {
                                        delx = nbrsLocMask[iT * 3]; dely = nbrsLocMask[iT * 3 + 1]; delz = nbrsLocMask[iT * 3 + 2];
                                        plint id1 = util::roundToInt(lattices[distLloc]->get(iXd + delx, iYd + dely, iZd + delz).computeDensity());
                                        if (id0 >= id1) { tmp1Vector.push_back(delx); tmp1Vector.push_back(dely); tmp1Vector.push_back(delz); ++tmp1Len; }
                                    }
                                    if (tmp1Len > 1) {
                                        plint randLoc = std::rand() % tmp1Len;
                                        delx = tmp1Vector[randLoc * 3]; dely = tmp1Vector[randLoc * 3 + 1]; delz = tmp1Vector[randLoc * 3 + 2];
                                    }
                                    else if (tmp1Len == 1) {
                                        delx = tmp1Vector[0]; dely = tmp1Vector[1]; delz = tmp1Vector[2];
                                    }
                                    else { // tmp1Len == 0, purely random redistribution
                                        plint randLoc = 0;
                                        if (nbrTlen > 1) { randLoc = std::rand() % nbrTlen; }
                                        delx = nbrsLocMask[randLoc * 3]; dely = nbrsLocMask[randLoc * 3 + 1]; delz = nbrsLocMask[randLoc * 3 + 2];

                                        // first check to see if bdry_dir is greater than 0
                                        if (bdry_dir > 0) {
                                          // of delx is 1 and bdry_dir is one of (1, 7, or 10), set push_dir to 1
                                          if (delx == 1 && (bdry_dir == 1 || bdry_dir == 7 || bdry_dir == 10)) { 
                                            push_dir = 1; 
                                          } 
                                          // if dely is 1 and bdry_dir is one of (2, 7, 8), set push_dir to 2
                                          else if (dely == 1 && (bdry_dir == 2 || bdry_dir == 7 || bdry_dir == 8)) { 
                                            push_dir = 2; 
                                          } 
                                          // if delz is 1 and bdry_dir is one of (3, 7, 8) set push_dir to 3
                                          else if (delz == 1 && (bdry_dir == 3 || bdry_dir == 7 || bdry_dir == 8)) { 
                                            push_dir = 3; 
                                          } 
                                          // if delx is -1 and bdry_dir is one of (4, 8, 9), set push_dir to 4
                                          else if (delx == -1 && (bdry_dir == 4 || bdry_dir == 8 || bdry_dir == 9)) { 
                                            push_dir = 4; 
                                          } 
                                          // if dely is -1 and bdry_dir is one of (5, 9, 10), set push_dir to 5
                                          else if (dely == -1 && (bdry_dir == 5 || bdry_dir == 9 || bdry_dir == 10)) {
                                            push_dir = 5; 
                                          } 
                                          // if delz is -1 and bdry_dir is one of (6, 9, 10) set push_dir to 6
                                          else if (delz == -1 && (bdry_dir == 6 || bdry_dir == 9 || bdry_dir == 10)) { 
                                            push_dir = 6; 
                                          } 
                                          // if none of the conditions are true, set push_dir to 0
                                          else { 
                                            push_dir = 0; 
                                          }
                                        }
                                        

                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            Array<T, 7> g;
                                            plint iXb1 = iX0 + vec_offset[iM].x, iYb1 = iY0 + vec_offset[iM].y, iZb1 = iZ0 + vec_offset[iM].z;; // original biomass lattice
                                            T shove_bmass = (lattices[iM]->get(iXb1, iYb1, iZb1).computeDensity()) / bMt * bMd;
                                            lattices[iM]->get(iXb1, iYb1, iZb1).getPopulations(g);
                                            g[0] -= (T)(shove_bmass) / 4; g[1] -= (T)(shove_bmass) / 8; g[2] -= (T)(shove_bmass) / 8; g[3] -= (T)(shove_bmass) / 8; g[4] -= (T)(shove_bmass) / 8; g[5] -= (T)(shove_bmass) / 8; g[6] -= (T)(shove_bmass) / 8;
                                            lattices[iM]->get(iXb1, iYb1, iZb1).setPopulations(g);
                                            // Update storage lattices if the selected neighbor location is outside of subdomain boundaries
                                            if (push_dir > 0) {
                                                plint iXb2 = iX0 + vec_offset[iM + numbM].x; plint iYb2 = iY0 + vec_offset[iM + numbM].y; plint iZb2 = iZ0 + vec_offset[iM + numbM].z;
                                                g[0] = (T)(shove_bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(shove_bmass - 1) / 8;
                                                lattices[iM + numbM]->get(iXb2, iYb2, iZb2).setPopulations(g);
                                                lattices[iM + numbM]->get(iXb2, iYb2, iZb2).getDynamics().setOmega(push_dir);
                                            }
                                            // Directly update biomass lattices if the selected neighbor location is inside of subdomain boundaries
                                            else {
                                                lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).getPopulations(g);
                                                g[0] += (T)(shove_bmass) / 4; g[1] += (T)(shove_bmass) / 8; g[2] += (T)(shove_bmass) / 8; g[3] += (T)(shove_bmass) / 8; g[4] += (T)(shove_bmass) / 8; g[5] += (T)(shove_bmass) / 8; g[6] += (T)(shove_bmass) / 8;
                                                lattices[iM]->get(iXb1 + delx, iYb1 + dely, iZb1 + delz).setPopulations(g);
                                            }
                                        }
                                    }
                                }


                            }

                        }
                    }
                }
            }
        }
            virtual BlockDomain::DomainT appliesTo() const {
                // Don't apply to envelope, because nearest neighbors need to be accessed.
                return BlockDomain::bulk;
            }
            virtual halfPushExcessBiomass3D<T, Descriptor>* clone() const {
                return new halfPushExcessBiomass3D<T, Descriptor>(*this);
            }
            void getTypeOfModification(std::vector<modif::ModifT>&modified) const {
                plint numbM = (length - 3) / 2;
                for (plint iB = 0; iB < numbM; ++iB) {
                    modified[iB] = modif::staticVariables;
                    modified[iB + numbM] = modif::allVariables;
                }
                modified[length - 3] = modif::nothing;
                modified[length - 2] = modif::nothing;
                modified[length - 1] = modif::nothing;
            }
    private:
        T Bmax;
        plint nx, ny, nz, bdryGap, length, solid, bb;
        std::vector<plint> pore;
     }; 



        // redistribute the excessive biomass (pull)
        template<typename T, template<typename U> class Descriptor>
        class pullExcessBiomass3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
        {
        public:
            pullExcessBiomass3D(plint nx_, plint ny_, plint nz_, plint bdryGap_, plint length_)
                : nx(nx_), ny(ny_), nz(nz_), bdryGap(bdryGap_), length(length_)
            {}
            // lattices[0~(#ofbM-1)] = original biomass lattices
            // lattices[#ofbM~(len-3)] = copy biomass lattices
            // lattices[len-3] = total biomass lattice
            // lattices[len-2] = mask lattice
            // lattices[len-1] = dist lattice
            virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
                std::vector<Dot3D> vec_offset;
                plint numbM = (length - 2) / 2;
                for (plint iL = 0; iL < length; ++iL) {
                    vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL]));
                }
                Dot3D absoluteOffset = lattices[0]->getLocation();
                plint iX0 = domain.x0, absX = iX0 + absoluteOffset.x;
                if (absX >= bdryGap) {
                    for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                        for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                            for (plint iM = 0; iM < numbM; ++iM) {
                                plint iXc = iX0 + vec_offset[iM + numbM].x, iYc = iY0 + vec_offset[iM + numbM].y, iZc = iZ0 + vec_offset[iM + numbM].z;
                                plint dir_id = util::roundToInt(lattices[iM + numbM]->get((iXc - 1), iYc, iZc).getDynamics().getOmega());
                                if (dir_id == 1) {
                                    T nbrbM = lattices[iM + numbM]->get((iXc - 1), iYc, iZc).computeDensity();
                                    if (nbrbM > thrd) {
                                        Array<T, 7> g;
                                        plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                        lattices[iM]->get(iXb, iYb, iZb).getPopulations(g);
                                        g[0] += (T)(nbrbM) / 4; g[1] += (T)(nbrbM) / 8; g[2] += (T)(nbrbM) / 8; g[3] += (T)(nbrbM) / 8; g[4] += (T)(nbrbM) / 8; g[5] += (T)(nbrbM) / 8; g[6] += (T)(nbrbM) / 8;
                                        lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                    }
                                }
                            }
                        }
                    }
                }
                iX0 = domain.x1; absX = iX0 + absoluteOffset.x;
                if (absX <= (nx - 1 - bdryGap)) {
                    for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                        for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                            for (plint iM = 0; iM < numbM; ++iM) {
                                plint iXc = iX0 + vec_offset[iM + numbM].x, iYc = iY0 + vec_offset[iM + numbM].y, iZc = iZ0 + vec_offset[iM + numbM].z;
                                plint dir_id = util::roundToInt(lattices[iM + numbM]->get((iXc + 1), iYc, iZc).getDynamics().getOmega());
                                if (dir_id == 3) {
                                    T nbrbM = lattices[iM + numbM]->get((iXc + 1), iYc, iZc).computeDensity();
                                    if (nbrbM > thrd) {
                                        Array<T, 7> g;
                                        plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                        lattices[iM]->get(iXb, iYb, iZb).getPopulations(g);
                                        g[0] += (T)(nbrbM) / 4; g[1] += (T)(nbrbM) / 8; g[2] += (T)(nbrbM) / 8; g[3] += (T)(nbrbM) / 8; g[4] += (T)(nbrbM) / 8; g[5] += (T)(nbrbM) / 8; g[6] += (T)(nbrbM) / 8;
                                        lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                    }
                                }
                            }
                        }
                    }
                }
                plint iY0 = domain.y0, absY = iY0 + absoluteOffset.y;
                if (absY > 0) {
                    for (iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                        for (plint iZ0 = domain.x0; iZ0 <= domain.z1; ++iZ0) {
                            for (plint iM = 0; iM < numbM; ++iM) {
                                plint iXc = iX0 + vec_offset[iM + numbM].x, iYc = iY0 + vec_offset[iM + numbM].y, iZc = iZ0 + vec_offset[iM + numbM].z;
                                plint dir_id = util::roundToInt(lattices[iM + numbM]->get(iXc, (iYc - 1), iZc).getDynamics().getOmega());
                                if (dir_id == 2) {
                                    T nbrbM = lattices[iM + numbM]->get(iXc, (iYc - 1), iZc).computeDensity();
                                    if (nbrbM > thrd) {
                                        Array<T, 7> g;
                                        plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                        lattices[iM]->get(iXb, iYb,iZb).getPopulations(g);
                                        g[0] += (T)(nbrbM) / 4; g[1] += (T)(nbrbM) / 8; g[2] += (T)(nbrbM) / 8; g[3] += (T)(nbrbM) / 8; g[4] += (T)(nbrbM) / 8; g[5] += (T)(nbrbM) / 8; g[6] += (T)(nbrbM) / 8;
                                        lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                    }
                                }
                            }
                        }
                    }
                }

                iY0 = domain.y1, absY = iY0 + absoluteOffset.y;
                if (absY < (ny - 1)) {
                    for (iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                        for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                            for (plint iM = 0; iM < numbM; ++iM) {
                                plint iXc = iX0 + vec_offset[iM + numbM].x, iYc = iY0 + vec_offset[iM + numbM].y, iZc = iZ0 + vec_offset[iM + numbM].z;
                                plint dir_id = util::roundToInt(lattices[iM + numbM]->get(iXc, (iYc + 1), iZc).getDynamics().getOmega());
                                if (dir_id == 4) {
                                    T nbrbM = lattices[iM + numbM]->get(iXc, (iYc + 1), iZc).computeDensity();
                                    if (nbrbM > thrd) {
                                        Array<T, 7> g;
                                        plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                        lattices[iM]->get(iXb, iYb, iZb).getPopulations(g);
                                        g[0] += (T)(nbrbM) / 4; g[1] += (T)(nbrbM) / 8; g[2] += (T)(nbrbM) / 8; g[3] += (T)(nbrbM) / 8; g[4] += (T)(nbrbM) / 8; g[5] += (T)(nbrbM) / 8; g[6] += (T)(nbrbM) / 8;
                                        lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                    }
                                }
                            }
                        }
                    }
                }

                plint iZ0 = domain.z0, absZ = iZ0 + absoluteOffset.z;
                if (absZ > 0) {
                    for (iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                        for (iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                            for (plint iM = 0; iM < numbM; ++iM) {
                                plint iXc = iX0 + vec_offset[iM + numbM].x, iYc = iY0 + vec_offset[iM + numbM].y, iZc = iZ0 + vec_offset[iM + numbM].z;
                                plint dir_id = util::roundToInt(lattices[iM + numbM]->get(iXc, iYc, (iZc - 1)).getDynamics().getOmega());
                                if (dir_id == 10) {
                                    T nbrbM = lattices[iM + numbM]->get(iXc, iYc, (iZc + 1)).computeDensity();
                                    if (nbrbM > thrd) {
                                        Array<T, 7> g;
                                        plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                        lattices[iM]->get(iXb, iYb, iZb).getPopulations(g);
                                        g[0] += (T)(nbrbM) / 4; g[1] += (T)(nbrbM) / 8; g[2] += (T)(nbrbM) / 8; g[3] += (T)(nbrbM) / 8; g[4] += (T)(nbrbM) / 8; g[5] += (T)(nbrbM) / 8; g[6] += (T)(nbrbM) / 8;
                                        lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                    }
                                }

                            }
                        }
                    }
                }

                iZ0 = domain.z1, absZ = iZ0 + absoluteOffset.z;
                if (absZ < (nz - 1)) {
                    for (iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                        for (iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                            for (plint iM = 0; iM < numbM; ++iM) {
                                plint iXc = iX0 + vec_offset[iM + numbM].x, iYc = iY0 + vec_offset[iM + numbM].y, iZc = iZ0 + vec_offset[iM + numbM].y;
                                plint dir_id = util::roundToInt(lattices[iM + numbM]->get(iXc, iYc, (iZc + 1)).getDynamics().getOmega());
                                if (dir_id == 9) {
                                    T nbrbM = lattices[iM + numbM]->get(iXc, iYc, (iZc - 1)).computeDensity();
                                    if (nbrbM > thrd) {
                                        Array<T, 7> g;
                                        plint iXb = iX0 + vec_offset[iM].x, iYb = iY0 + vec_offset[iM].y, iZb = iZ0 + vec_offset[iM].z;
                                        lattices[iM]->get(iXb, iYb, iZb).getPopulations(g);
                                        g[0] += (T)(nbrbM) / 4; g[1] += (T)(nbrbM) / 8; g[2] += (T)(nbrbM) / 8; g[3] += (T)(nbrbM) / 8; g[4] += (T)(nbrbM) / 8; g[5] += (T)(nbrbM) / 8; g[6] += (T)(nbrbM) / 8;
                                        lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                    }
                                }

                            }
                        }
                    }
                }
            }
            virtual BlockDomain::DomainT appliesTo() const {
                // Don't apply to envelope, because nearest neighbors need to be accessed.
                return BlockDomain::bulk;
            }
            virtual pullExcessBiomass3D<T, Descriptor>* clone() const {
                return new pullExcessBiomass3D<T, Descriptor>(*this);
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                plint numbM = (length - 3) / 2;
                for (plint iB = 0; iB < numbM; ++iB) {
                    modified[iB] = modif::staticVariables;
                    modified[iB + numbM] = modif::nothing;
                }
                modified[length - 3] = modif::nothing;
                modified[length - 2] = modif::nothing;
                modified[length - 1] = modif::nothing;
            }
        private:
            plint nx, ny, nz, bdryGap, length;
        };
        

        // Update nonlocal mask number at every timestep
        template<typename T, template<typename U> class Descriptor>
        class updateLocalMaskNtotalLattices3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
        {
        public:
            updateLocalMaskNtotalLattices3D(plint nx_, plint ny_, plint nz_, plint length_, plint bb_, plint solid_, std::vector<std::vector<plint>> bio_, std::vector<plint> pore_, T bMassFrac_, T maxbMassRho_)
                : nx(nx_), ny(ny_), nz(nz_), length(length_), bb(bb_), solid(solid_), bio(bio_), pore(pore_), thrdbMassRho(bMassFrac_* maxbMassRho_)
            {}
            // lattices[0~(#ofbM-1)] = original biomass lattices
            // lattices[#ofbM~(len-3)] = copy biomass lattices
            // lattices[len-3] = total biomass lattice
            // lattices[len-2] = mask lattice
            // lattices[len-1] = age lattice
            virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
                plint maskLloc = length - 2, bMtLloc = length - 3;
                std::vector<Dot3D> vec_offset;
                Dot3D absoluteOffset = lattices[0]->getLocation();
                for (plint iL = 0; iL < length; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL])); }
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iXm = iX + vec_offset[maskLloc].x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iYm = iY + vec_offset[maskLloc].y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZm = iZ + vec_offset[maskLloc].z;
                            plint mask = util::roundToInt(lattices[maskLloc]->get(iXm, iYm, iZm).computeDensity());
                            if (mask != bb && mask != solid) {
                                T bmass = 0.; plint newmask = 0;
                                T bMt = lattices[bMtLloc]->get(iX + vec_offset[bMtLloc].x, iY + vec_offset[bMtLloc].y, iZ + vec_offset[bMtLloc].z).computeDensity();
                                for (size_t iM = 0; iM < bio.size(); ++iM) {
                                    T bm = lattices[iM]->get(iX + vec_offset[iM].x, iY + vec_offset[iM].y, iZ + vec_offset[iM].z).computeDensity();
                                    if (bm > thrd) { bmass += bm; newmask += bio[iM][0]; }
                                }
                                // update total biomass density
                                if (std::abs(bmass - bMt) > thrd) {
                                    Array<T, 7> g;
                                    g[0] = (T)(bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(bmass - 1) / 8;
                                    plint iXt = iX + vec_offset[bMtLloc].x, iYt = iY + vec_offset[bMtLloc].y, iZt = iZ + vec_offset[bMtLloc].z;
                                    lattices[bMtLloc]->get(iXt, iYt, iZt).setPopulations(g);
                                }
                                // update mask number
                                bool poreflag = 0;
                                for (size_t iP = 0; iP < pore.size(); ++iP) {
                                    if (mask == pore[iP]) { poreflag = 1; break; }
                                }
                                if (poreflag == 0 && bmass < thrdbMassRho) {
                                    newmask = pore[0];
                                    if (pore.size() > 1) {
                                        plint absX = iX + absoluteOffset.x, absY = iY + absoluteOffset.y, absZ = iZ + absoluteOffset.z, nbrs = 8;
                                        std::vector<plint> delXYZ;
                                        if (absX == 0 && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {
                                            nbrs = 3; // nbrs is the number of neighbors depending on the current location
                                            delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                        }
                                        else if (absX == (nx - 1) && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {
                                            nbrs = 3;
                                            delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);
                                            
                                        }
                                        else if (absX > 0 && absX < (nx - 1) && absY == 0 && absZ == 0) {
                                            nbrs = 3;
                                            delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);
                                            
                                        }
                                        else if (absX > 0 && absX < (nx - 1) && absY == (ny - 1) && absZ == (nz - 1)) {
                                            nbrs = 3;
                                            delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                            
                                        }
                                        else if (absX == 0 && absY == 0 && absZ == 0) {
                                            nbrs = 2;
                                            delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);
                                        }
                                            
                                        else if ( absX == 0 && absY == (ny - 1) && absZ == (nz - 1)) {
                                            nbrs = 2;
                                            delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                           
                                        }
                                        else if (absX == (nx - 1) && absY == 0 && absZ == 0) {
                                            nbrs = 2;
                                            delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);
                                           
                                        }
                                        else if (absX == (nx - 1) && absY == ny - 1 && absZ == nz - 1) {
                                            nbrs = 2;
                                            delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                        }
                                        else {
                                            nbrs = 4;
                                            delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                            delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);
                                        }
                                        for (plint iT = 0; iT < nbrs; ++iT) {
                                            plint delx = delXYZ[iT * 2], dely = delXYZ[iT * 2 + 1], delz = delXYZ[iT * 2 + 1];
                                            plint nbrmask = util::roundToInt(lattices[maskLloc]->get(iXm + delx, iYm + dely, iZm + delz).computeDensity());
                                            if (nbrmask == bb) {
                                                newmask = pore[1];
                                            }
                                        }
                                    }
                                    Array<T, 7> g;
                                    g[0] = (T)(newmask - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(newmask - 1) / 8;
                                    lattices[maskLloc]->get(iXm, iYm, iZm).setPopulations(g);
                                }
                                else if (poreflag == 1 && bmass >= thrdbMassRho) {
                                    if (newmask > 0) {
                                        Array<T, 7> g;
                                        g[0] = (T)(newmask - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(newmask - 1) / 8;
                                        lattices[maskLloc]->get(iXm, iYm, iZm).setPopulations(g);
                                    }
                                    else {
                                        std::cout << "Error: Updating mask failed.\n";
                                        exit(EXIT_FAILURE);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulk;
            }
            virtual updateLocalMaskNtotalLattices3D<T, Descriptor>* clone() const {
                return new updateLocalMaskNtotalLattices3D<T, Descriptor>(*this);
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                for (plint iL = 0; iL < length - 3; ++iL) {
                    modified[iL] = modif::nothing;
                }
                modified[length - 3] = modif::staticVariables;
                modified[length - 2] = modif::staticVariables;
                modified[length - 1] = modif::nothing;
            }
        private:
            plint nx, ny, nz, length, bb, solid;
            std::vector<std::vector<plint>> bio;
            std::vector<plint> pore;
            T thrdbMassRho;
        };

        

        template<typename T, template<typename U> class Descriptor>
        class fdDiffusion3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
        {
        public:
            fdDiffusion3D(plint nx_, plint ny_, plint nz_, plint length_, plint bdryGap_, T nu_)
                : nx(nx_), ny(ny_), nz(nz_), length(length_), bdryGap(bdryGap_), nu(nu_)
            {}
            // lattices[0~(#ofbM-1)] = original biomass lattices
            // lattices[#ofbM~(len-2)] = copy biomass lattices
            // lattices[len-1] = mask lattice
            virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
                std::vector<Dot3D> vec_offset;
                plint maskLloc = length - 1, numbM = (length - 1) / 2;
                // plint bMtLloc = length-2;
                for (plint iL = 0; iL < length; ++iL) {
                    vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL]));
                }
                Dot3D absoluteOffset = lattices[0]->getLocation();
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint absX = iX + absoluteOffset.x;
                    if (absX >= bdryGap && absX <= nx - 1 - bdryGap) {
                        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                            plint absY = iY + absoluteOffset.y;
                            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                                plint absZ = iZ + absoluteOffset.z;
                                // collect mask numbers
                                plint iXm = iX + vec_offset[maskLloc].x, iYm = iY + vec_offset[maskLloc].y, iZm = iZ + vec_offset[maskLloc].z;
                                plint mask = util::roundToInt(lattices[maskLloc]->get(iXm, iYm, iZm).computeDensity());
                                if (mask > 1) {
                                    plint mxp = 0, mxn = 0, myp = 0, myn = 0, mzp = 0, mzn = 0;
                                    if (absX == bdryGap) { mxp = lattices[maskLloc]->get(iXm + 1, iYm, iZm).computeDensity(); }
                                    else if (absX == nx - 1 - bdryGap) { mxn = lattices[maskLloc]->get(iXm - 1, iYm, iZm).computeDensity(); }
                                    else { mxp = lattices[maskLloc]->get(iXm + 1, iYm, iZm).computeDensity(); mxn = lattices[maskLloc]->get(iXm - 1, iYm, iZm).computeDensity(); }

                                    if (absY == 0) { myp = lattices[maskLloc]->get(iXm, iYm + 1,iZm ).computeDensity(); }
                                    else if (absY == ny - 1) { myn = lattices[maskLloc]->get(iXm, iYm - 1, iZm).computeDensity(); }
                                    else { myp = lattices[maskLloc]->get(iXm, iYm + 1, iZm).computeDensity(); myn = lattices[maskLloc]->get(iXm, iYm - 1, iZm).computeDensity(); }

                                    if (absZ == 0) { mzp = lattices[maskLloc]->get(iXm, iYm, iZm + 1).computeDensity(); }
                                    else if (absZ == nz - 1) { mzn = lattices[maskLloc]->get(iXm, iYm, iZm - 1).computeDensity(); }
                                    else { mzp = lattices[maskLloc]->get(iXm, iYm, iZm + 1).computeDensity(); mzn = lattices[maskLloc]->get(iXm, iYm, iZm - 1).computeDensity(); }

                                    // collect total biomass densities
                                    // plint iXt = iX + vec_offset[bMtLloc].x, iYt = iY + vec_offset[bMtLloc].y;
                                    // T bMt0 = lattices[bMtLloc]->get(iXt,iYt).computeDensity();
                                    // T bMtxp=0,bMtxn=0,bMtyp=0,bMtyn=0;
                                    std::vector<T> bxp(numbM, 0), bxn(numbM, 0), byp(numbM, 0), byn(numbM, 0), bzp(numbM, 0), bzn(numbM, 0), b0(numbM, 0);
                                    for (plint iM = 0; iM < numbM; ++iM) {
                                        plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                        b0[iM] = lattices[iM + numbM]->get(iXb, iYb, iZb).computeDensity();
                                    }
                                    if (absX == bdryGap || mxn < 2) {
                                        // bMtxp=lattices[bMtLloc]->get(iXt+1,iYt).computeDensity();
                                        // bMtxn=bMt0;
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            bxp[iM] = lattices[iM + numbM]->get(iXb + 1, iYb, iZb).computeDensity();
                                            bxn[iM] = b0[iM];
                                        }
                                    }
                                    else if (absX == nx - 1 - bdryGap || mxp < 2) {
                                        // bMtxp=bMt0;
                                        // bMtxn=lattices[bMtLloc]->get(iXt-1,iYt).computeDensity();
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            bxp[iM] = b0[iM];
                                            bxn[iM] = lattices[iM + numbM]->get(iXb - 1, iYb, iZb).computeDensity();
                                        }
                                    }
                                    else {
                                        // bMtxp=lattices[bMtLloc]->get(iXt+1,iYt).computeDensity();
                                        // bMtxn=lattices[bMtLloc]->get(iXt-1,iYt).computeDensity();
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            bxp[iM] = lattices[iM + numbM]->get(iXb + 1, iYb, iZb).computeDensity();
                                            bxn[iM] = lattices[iM + numbM]->get(iXb - 1, iYb, iZb).computeDensity();
                                        }
                                    }
                                    if (absY == 0 || myn < 2) {
                                        // bMtyp=lattices[bMtLloc]->get(iXt,iYt+1).computeDensity();
                                        // bMtyn=bMt0;
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            byp[iM] = lattices[iM + numbM]->get(iXb, iYb + 1, iZb).computeDensity();
                                            byn[iM] = b0[iM];
                                        }
                                    }
                                    else if (absY == ny - 1 || myp < 2) {
                                        // bMtyp=bMt0;
                                        // bMtyn=lattices[bMtLloc]->get(iXt,iYt-1).computeDensity();
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            byp[iM] = b0[iM];
                                            byn[iM] = lattices[iM + numbM]->get(iXb, iYb - 1, iZb).computeDensity();
                                        }
                                    }
                                    else {
                                        // bMtyp=lattices[bMtLloc]->get(iXt,iYt+1).computeDensity();
                                        // bMtyn=lattices[bMtLloc]->get(iXt,iYt-1).computeDensity();
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            byp[iM] = lattices[iM + numbM]->get(iXb, iYb + 1, iZb).computeDensity();
                                            byn[iM] = lattices[iM + numbM]->get(iXb, iYb - 1, iZb).computeDensity();
                                        }
                                    }
                                    if (absZ == 0 || mzn < 2) {
                                        // bMtyp=lattices[bMtLloc]->get(iXt,iYt+1).computeDensity();
                                        // bMtyn=bMt0;
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            bzp[iM] = lattices[iM + numbM]->get(iXb, iYb, iZb + 1).computeDensity();
                                            bzn[iM] = b0[iM];
                                        }
                                    }
                                    else if (absZ == nz - 1 || mzp < 2) {
                                        // bMtyp=bMt0;
                                        // bMtyn=lattices[bMtLloc]->get(iXt,iYt-1).computeDensity();
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            bzp[iM] = b0[iM];
                                            bzn[iM] = lattices[iM + numbM]->get(iXb, iYb, iZb - 1).computeDensity();
                                        }
                                    }
                                    else {
                                        // bMtyp=lattices[bMtLloc]->get(iXt,iYt+1).computeDensity();
                                        // bMtyn=lattices[bMtLloc]->get(iXt,iYt-1).computeDensity();
                                        for (plint iM = 0; iM < numbM; ++iM) {
                                            plint iXb = iX + vec_offset[iM + numbM].x, iYb = iY + vec_offset[iM + numbM].y, iZb = iZ + vec_offset[iM + numbM].z;
                                            bzp[iM] = lattices[iM + numbM]->get(iXb, iYb, iZb + 1).computeDensity();
                                            bzn[iM] = lattices[iM + numbM]->get(iXb, iYb, iZb + 1).computeDensity();
                                        }
                                    }
                                    // if (bMt0 < thrd) bMt0=0;
                                    // if (bMtxp < thrd) bMtxp=0;
                                    // if (bMtxn < thrd) bMtxn=0;
                                    // if (bMtyp < thrd) bMtyp=0;
                                    // if (bMtyn < thrd) bMtyn=0;
                                    for (plint iM = 0; iM < numbM; ++iM) {
                                        if (b0[iM] < thrd) b0[iM] = 0;
                                        if (bxp[iM] < thrd) bxp[iM] = 0;
                                        if (bxn[iM] < thrd) bxn[iM] = 0;
                                        if (byp[iM] < thrd) byp[iM] = 0;
                                        if (byn[iM] < thrd) byn[iM] = 0;
                                        if (byp[iM] < thrd) bzp[iM] = 0;
                                        if (byn[iM] < thrd) bzn[iM] = 0;
                                        T new_bM = b0[iM] + nu * ((bxp[iM] - 2 * b0[iM] + bxn[iM]) + (byp[iM] - 2 * b0[iM] + byn[iM]) + (bzp[iM] - 2 * b0[iM] + bzn[iM]));
                                        if (new_bM > thrd) {
                                            Array<T, 7> g;
                                            plint iXb = iX + vec_offset[iM].x, iYb = iY + vec_offset[iM].y, iZb = iZ + vec_offset[iM].z;
                                            g[0] = (T)(new_bM - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(new_bM - 1) / 8;
                                            lattices[iM]->get(iXb, iYb, iZb).setPopulations(g);
                                        }
                                    }

                                    // T bMt1 = bMt0 + nu*( (bMtxp-2*bMt0+bMtxn) + (bMtyp-2*bMt0+bMtyn) );
                                    // T delta_bMt = bMt1 - bMt0;
                                    // if (delta_bMt > thrd) {
                                    //     for (plint iB=0; iB<numbM; ++iB) {
                                    //         plint iXb0 = iX + vec_offset[iB+numbM].x, iYb0 = iY + vec_offset[iB+numbM].y;
                                    //         T old_bM = lattices[iB+numbM]->get(iXb0,iYb0).computeDensity();
                                    //         T ratio_bM = delta_bMt*old_bM/bMt0;
                                    //         Array<T,5> g;
                                    //         plint iXb1 = iX + vec_offset[iB].x, iYb1 = iY + vec_offset[iB].y;
                                    //         lattices[iB+numbM]->get(iXb1,iYb1).getPopulations(g);
                                    //         g[0]+=(T) (ratio_bM)/3; g[1]+=(T) (ratio_bM)/6; g[2]+=(T) (ratio_bM)/6; g[3]+=(T) (ratio_bM)/6; g[4]+=(T) (ratio_bM)/6;
                                    //         lattices[iB]->get(iXb1,iYb1).setPopulations(g);
                                    //     }
                                    // }
                                }
                            }
                        }
                    }
                }
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulk;
            }
            virtual fdDiffusion3D<T, Descriptor>* clone() const {
                return new fdDiffusion3D<T, Descriptor>(*this);
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                plint numbM = (length - 1) / 2;
                for (plint iB = 0; iB < numbM; ++iB) {
                    modified[iB] = modif::staticVariables;
                    modified[iB + numbM] = modif::nothing;
                }
                modified[length - 1] = modif::nothing;
            }
        private:
            plint nx, ny, nz, length, bdryGap;
            T nu;
        };
        

        // update solute diffusivity at a biomass voxel
        template<typename T, template<typename U> class Descriptor>
        class updateSoluteDynamics3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
        {
        public:
            updateSoluteDynamics3D(plint subsNum_, plint bb_, plint solid_, std::vector<plint> pore_, std::vector<T> substrOMEGAinbMass_, std::vector<T> substrOMEGAinPore_)
                : subsNum(subsNum_), bb(bb_), solid(solid_), pore(pore_), substrOMEGAinbMass(substrOMEGAinbMass_), substrOMEGAinPore(substrOMEGAinPore_)
            {}
            // lattices[0~(#ofSubs-1)] = substrate lattices
            // lattices[#ofSubs] = mask lattice
            virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
                std::vector<Dot3D> vec_offset;
                for (plint iT = 0; iT < subsNum + 1; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT])); }
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iXm = iX + vec_offset[subsNum].x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iYm = iY + vec_offset[subsNum].y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZm = iZ + vec_offset[subsNum].z;
                            plint mask = util::roundToInt(lattices[subsNum]->get(iXm, iYm, iZm).computeDensity());
                            if (mask != bb && mask != solid) {
                                bool poreflag = 0;
                                for (size_t iP = 0; iP < pore.size(); ++iP) { if (mask == pore[iP]) { poreflag = 1; break; } }
                                for (plint iS = 0; iS < subsNum; ++iS) {
                                    plint iXs = iX + vec_offset[iS].x, iYs = iY + vec_offset[iS].y, iZs = iZ + vec_offset[iS].z;
                                    T omega = lattices[iS]->get(iXs, iYs, iZs).getDynamics().getOmega();
                                    T bMassOmega = substrOMEGAinbMass[iS];
                                    if (poreflag == 0 && std::abs(omega - bMassOmega) > thrd) {
                                        lattices[iS]->get(iXs, iYs, iZs).getDynamics().setOmega(bMassOmega); // pore to bfilm
                                    }
                                    else if (poreflag == 1 && std::abs(omega - bMassOmega) < thrd) {
                                        lattices[iS]->get(iXs, iYs, iZs).getDynamics().setOmega(substrOMEGAinPore[iS]); // bfilm to pore
                                    }
                                }
                            }
                        }
                    }
                }
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            virtual updateSoluteDynamics3D<T, Descriptor>* clone() const {
                return new updateSoluteDynamics3D<T, Descriptor>(*this);
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                for (plint iS = 0; iS < subsNum; ++iS) {
                    modified[iS] = modif::dynamicVariables;
                }
                modified[subsNum] = modif::nothing;
            }
        private:
            plint subsNum, bb, solid;
            std::vector<plint> pore;
            std::vector<T> substrOMEGAinbMass, substrOMEGAinPore;
        };







        // update planktonic biomass diffusivity at a biomass voxel
        template<typename T, template<typename U> class Descriptor>
        class updateBiomassDynamics3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
        {
        public:
            updateBiomassDynamics3D(plint bioNum_, plint bb_, plint solid_, std::vector<plint> pore_, std::vector<T> bmassOMEGAinbMass_, std::vector<T> bmassOMEGAinPore_)
                : bioNum(bioNum_), bb(bb_), solid(solid_), pore(pore_), bmassOMEGAinbMass(bmassOMEGAinbMass_), bmassOMEGAinPore(bmassOMEGAinPore_)
            {}
            // lattices[0~(#ofbMs-1)] = planktonic biomass lattices
            // lattices[#ofbMs] = mask lattice
            virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
                std::vector<Dot3D> vec_offset;
                for (plint iT = 0; iT < bioNum + 1; ++iT) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iT])); }
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iXm = iX + vec_offset[bioNum].x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iYm = iY + vec_offset[bioNum].y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZm = iZ + vec_offset[bioNum].z;
                            plint mask = util::roundToInt(lattices[bioNum]->get(iXm, iYm, iZm).computeDensity());
                            if (mask != bb && mask != solid) {
                                bool poreflag = 0;
                                for (size_t iP = 0; iP < pore.size(); ++iP) { if (mask == pore[iP]) { poreflag = 1; break; } }
                                for (plint iB = 0; iB < bioNum; ++iB) {
                                    plint iXb = iX + vec_offset[iB].x, iYb = iY + vec_offset[iB].y, iZb = iZ + vec_offset[iB].z;
                                    T omega = lattices[iB]->get(iXb, iYb, iZb).getDynamics().getOmega();
                                    T bMassOmega = bmassOMEGAinbMass[iB];
                                    if (poreflag == 0 && std::abs(omega - bMassOmega) > thrd) {
                                        lattices[iB]->get(iXb, iYb, iZb).getDynamics().setOmega(bMassOmega); // pore to bfilm
                                    }
                                    else if (poreflag == 1 && std::abs(omega - bMassOmega) < thrd) {
                                        lattices[iB]->get(iXb, iYb, iZb).getDynamics().setOmega(bmassOMEGAinPore[iB]); // bfilm to pore
                                    }
                                }
                            }
                        }
                    }
                }
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            virtual updateBiomassDynamics3D<T, Descriptor>* clone() const {
                return new updateBiomassDynamics3D<T, Descriptor>(*this);
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                for (plint iB = 0; iB < bioNum; ++iB) {
                    modified[iB] = modif::dynamicVariables;
                }
                modified[bioNum] = modif::nothing;
            }
        private:
            plint bioNum, bb, solid;
            std::vector<plint> pore;
            std::vector<T> bmassOMEGAinbMass, bmassOMEGAinPore;
        };

        // update the flow field
        template<typename T1, template<typename U> class Descriptor1, typename T2, template<typename U> class Descriptor2>
        class updateNsLatticesDynamics3D : public BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>
        {
        public:
            updateNsLatticesDynamics3D(T nsOmega_, T bioX_, std::vector<plint> pore_, plint solid_, plint bb_)
                : nsOmega(nsOmega_), bioX(bioX_), pore(pore_), solid(solid_), bb(bb_)
            {}
            // lattice0 = flow field lattice
            // lattice1 = mask field lattice
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice0, BlockLattice3D<T2, Descriptor2>& lattice1) {
                T bioOmega = 1 / (bioX * (1 / nsOmega - .5) + .5);
                Dot3D offset_12 = computeRelativeDisplacement(lattice0, lattice1);
                for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                    plint iX1 = iX0 + offset_12.x;
                    for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                        plint iY1 = iY0 + offset_12.y;
                        for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                            plint iZ1 = iZ0 + offset_12.z;
                            plint mask = util::roundToInt(lattice1.get(iX1, iY1, iZ1).computeDensity());
                            T currentOmega = lattice0.get(iX0, iY0, iZ0).getDynamics().getOmega();
                            if (mask != bb || mask != solid) {
                                bool poreflag = 0;
                                for (size_t iP = 0; iP < pore.size(); ++iP) { if (mask == pore[iP]) { poreflag = 1; break; } }
                                // from pore to biomass
                                if (poreflag == 0 && std::abs(currentOmega - nsOmega) < thrd) {
                                    if (bioX <= thrd) { lattice0.attributeDynamics(iX0, iY0, iZ0, new BounceBack<T1, Descriptor1>()); }
                                    else { lattice0.attributeDynamics(iX0, iY0, iZ0, new IncBGKdynamics<T1, Descriptor1>(bioOmega)); }
                                }
                                // from biomass to pore
                                else if (poreflag == 1 && std::abs(currentOmega - nsOmega) > thrd) {
                                    lattice0.attributeDynamics(iX0, iY0, iZ0, new IncBGKdynamics<T1, Descriptor1>(nsOmega));
                                    // lattice0.attributeDynamics(iX0, iY0, new IncMRTdynamics<T1,Descriptor1>(nsOmega) );
                                }
                            }
                        }
                    }
                }
            }
            virtual updateNsLatticesDynamics3D<T1, Descriptor1, T2, Descriptor2>* clone() const {
                return new updateNsLatticesDynamics3D<T1, Descriptor1, T2, Descriptor2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::dataStructure;
                modified[1] = modif::nothing;
            }
        private:
            T nsOmega, bioX;
            std::vector<plint> pore;
            plint solid, bb;
        };



  // Link geometry scalar numbers and maskLattice
        template<typename T1, template<typename U> class Descriptor, typename T2>
        class CopyGeometryScalar2maskLattice3D : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
        {
        public:
            CopyGeometryScalar2maskLattice3D(std::vector< std::vector<plint> > mask0_) : mask0(mask0_)
            {}
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
                Dot3D offset = computeRelativeDisplacement(lattice, field);
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iX1 = iX + offset.x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iY1 = iY + offset.y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZ1 = iZ + offset.z;
                            bool flag = 0; T mask1 = field.get(iX1, iY1, iZ1); T mask2 = -1;
                            for (size_t iM = 0; iM < mask0.size(); ++iM) {
                                for (size_t iN = 0; iN < mask0[iM].size(); ++iN) {
                                    if (mask1 == mask0[iM][iN]) {
                                        flag = 1;
                                        mask2 = mask0[iM][0];
                                        break;
                                    }
                                }
                                if (flag == 1) {
                                    break;
                                }
                            }
                            Array<T, 7> g;
                            if (flag == 0) {
                                g[0] = (T)(mask1 - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(mask1 - 1) / 8;
                            }
                            else {
                                g[0] = (T)(mask2 - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(mask2 - 1) / 8;
                            }
                            lattice.get(iX, iY, iZ).setPopulations(g); // allocate the mask number
                        }
                    }
                }
            }
            virtual CopyGeometryScalar2maskLattice3D<T1, Descriptor, T2>* clone() const {
                return new CopyGeometryScalar2maskLattice3D<T1, Descriptor, T2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::staticVariables;
                modified[1] = modif::nothing;
            }
        private:
            std::vector<std::vector<plint>> mask0;
        };

        // Link geometry scalar numbers and maskLattice
        template<typename T1, template<typename U> class Descriptor, typename T2>
        class CopyGeometryScalar2ageLattice3D : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
        {
        public:
            CopyGeometryScalar2ageLattice3D()
            {}
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) 
            {
                Dot3D offset = computeRelativeDisplacement(lattice, field);
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) 
                {
                    plint iX1 = iX + offset.x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) 
                    {
                        plint iY1 = iY + offset.y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) 
                        {
                            plint iZ1 = iZ + offset.z;
                            plint mask = field.get(iX1, iY1, iZ1);
                            if (mask < 0) { mask = -1; }
                            Array<T, 7> g;
                            g[0] = (T)(mask - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(mask - 1) / 8;
                            lattice.get(iX, iY, iZ).setPopulations(g); // allocate the mask number
                        }
                    }
                }
            }
            virtual CopyGeometryScalar2ageLattice3D<T1, Descriptor, T2>* clone() const {
                return new CopyGeometryScalar2ageLattice3D<T1, Descriptor, T2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::staticVariables;
                modified[1] = modif::nothing;
            }
        private:
        };

        // Link geometry scalar numbers and maskLattice
        template<typename T1, template<typename U> class Descriptor, typename T2>
        class CopyGeometryScalar2distLattice3D : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
        {
        public:
            CopyGeometryScalar2distLattice3D()
            {}
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
                Dot3D offset = computeRelativeDisplacement(lattice, field);
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iX1 = iX + offset.x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iY1 = iY + offset.y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZ1 = iZ + offset.z;
                            plint dist = field.get(iX1, iY1, iZ1);
                            Array<T, 7> g;
                            g[0] = (T)(dist - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(dist - 1) / 8;
                            lattice.get(iX, iY, iZ).setPopulations(g);
                        }
                    }
                }
            }
            virtual CopyGeometryScalar2distLattice3D<T1, Descriptor, T2>* clone() const {
                return new CopyGeometryScalar2distLattice3D<T1, Descriptor, T2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::staticVariables;
                modified[1] = modif::nothing;
            }
        private:
        };

        // Copy maskLattice to geometry field
        template<typename T1, template<typename U> class Descriptor, typename T2>
        class CopyLattice2ScalarField : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
        {
        public:
            CopyLattice2ScalarField()
            {}
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
                Dot3D offset = computeRelativeDisplacement(lattice, field);
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iX1 = iX + offset.x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iY1 = iY + offset.y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZ1 = iZ + offset.z;
                            field.get(iX1, iY1, iZ1) = util::roundToInt(lattice.get(iX, iY, iZ).computeDensity());
                        }
                    }
                }
            }
            virtual CopyLattice2ScalarField<T1, Descriptor, T2>* clone() const {
                return new CopyLattice2ScalarField<T1, Descriptor, T2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::nothing;
                modified[1] = modif::staticVariables;
            }
        private:
            std::vector<std::vector<plint>> mask0;
        };

        // initialize scalar biomass lattice
        template<typename T1, template<typename U> class Descriptor, typename T2>
        class initializeScalarLattice : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
        {
        public:
            initializeScalarLattice(std::vector<T> b0_, std::vector<plint> mask0_) : b0(b0_), mask0(mask0_)
            {}
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
                Dot3D offset = computeRelativeDisplacement(lattice, field);
                if (b0.size() != mask0.size()) {
                    std::cout << "ERROR: the size of vectors b0 and mask0 should be the same.\n";
                    exit(EXIT_FAILURE);
                }
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iX1 = iX + offset.x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iY1 = iY + offset.y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZ1 = iZ + offset.z;
                            T2 mask1 = field.get(iX1, iY1, iZ1);
                            for (size_t iN = 0; iN < b0.size(); ++iN) {
                                if (mask1 == mask0[iN]) {
                                    Array<T, 7> g;
                                    g[0] = (T)(b0[iN] - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(b0[iN] - 1) / 8;
                                    lattice.get(iX, iY, iZ).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
            virtual initializeScalarLattice<T1, Descriptor, T2>* clone() const {
                return new initializeScalarLattice<T1, Descriptor, T2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::staticVariables;
                modified[1] = modif::nothing;
            }
        private:
            std::vector<T> b0;
            std::vector<plint> mask0;
        };

        // initialize scalar biomass lattice
        template<typename T1, template<typename U> class Descriptor, typename T2>
        class stabilizeADElattice : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
        {
        public:
            stabilizeADElattice(T c0_, std::vector<plint> pore_, std::vector< std::vector<plint> > bio_) : c0(c0_), pore(pore_), bio(bio_)
            {}
            virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
                Dot3D offset = computeRelativeDisplacement(lattice, field);
                for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                    plint iX1 = iX + offset.x;
                    for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                        plint iY1 = iY + offset.y;
                        for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                            plint iZ1 = iZ + offset.z;
                            T2 mask = field.get(iX1, iY1, iZ1);
                            bool chk = 0;
                            for (size_t iP = 0; iP < pore.size(); ++iP) {
                                if (mask == pore[iP]) { chk = 1; break; }
                            }
                            for (size_t iB0 = 0; iB0 < bio.size(); ++iB0) {
                                if (chk == 1) { break; }
                                for (size_t iB1 = 0; iB1 < bio[iB0].size(); ++iB1) {
                                    if (mask == bio[iB0][iB1]) { chk = 1; break; }
                                }
                            }
                            if (chk == 1) {
                                if (c0<thrd && c0>-thrd) { c0 = 0; }
                                Array<T, 7> g;
                                g[0] = (T)(c0 - 1) / 4; g[1] = (T)(c0 - 1) / 8; g[2] = (T)(c0 - 1) / 8; g[3] = (T)(c0 - 1) / 8; g[4] = (T)(c0 - 1) / 8; g[5] = (T)(c0 - 1) / 8; g[6] = (T)(c0 - 1) / 8;
                                lattice.get(iX, iY, iZ).setPopulations(g);

                            }
                        }
                    }
                }
            }
            virtual stabilizeADElattice<T1, Descriptor, T2>* clone() const {
                return new stabilizeADElattice<T1, Descriptor, T2>(*this);
            }
            virtual BlockDomain::DomainT appliesTo() const {
                return BlockDomain::bulkAndEnvelope;
            }
            void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                modified[0] = modif::staticVariables;
                modified[1] = modif::nothing;
            }
        private:
            T c0;
            std::vector<plint> pore;
            std::vector< std::vector<plint> > bio;
        };


        // redefine age lattice
        template<typename T, template<typename U> class Descriptor>
        class updateAgeDistance3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
        {
        public:
            updateAgeDistance3D(T Bmax_, plint nx_, plint ny_, plint nz_)
                : Bmax(Bmax_), nx(nx_), ny(ny_), nz(nz_)
            {}
            // lattices[0] = age 0 lattices
            // lattices[1] = dist lattices
            // lattices[2] = total biomass lattice
            virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
                std::vector<Dot3D> vec_offset;
                Dot3D absoluteOffset = lattices[0]->getLocation();
                plint ageLloc = 0, distLloc = 1, bMtLloc = 2;
                for (plint iL = 0; iL < 3; ++iL) { vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL])); }
                for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
                    plint iXb = iX0 + vec_offset[bMtLloc].x;
                    for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                        plint iYb = iY0 + vec_offset[bMtLloc].y;
                        for (plint iZ0 = domain.y0; iZ0 <= domain.z1; ++iZ0) {
                            plint iZb = iZ0 + vec_offset[bMtLloc].z;
                            T B = lattices[bMtLloc]->get(iXb, iYb, iZb).computeDensity();
                            if (B > thrd) {
                                plint absX = iX0 + absoluteOffset.x, absY = iY0 + absoluteOffset.y, absZ = iZ0 + absoluteOffset.z;
                                plint iXa = iX0 + vec_offset[ageLloc].x, iYa = iY0 + vec_offset[ageLloc].y, iZa = iZ0 + vec_offset[ageLloc].z;
                                plint iXd = iX0 + vec_offset[distLloc].x, iYd = iY0 + vec_offset[distLloc].y, iZd = iZ0 + vec_offset[distLloc].z;
                                std::vector<plint> delXYZ; plint nbrs = 0;


                                if (absX == 0 && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {

                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == (nx - 1) && absY > 0 && absY < (ny - 1) && absZ > 0 && absZ < (nz - 1)) {

                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX > 0 && absX < (nx - 1) && absY == 0 && absZ == 0) {

                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX > 0 && absX < (nx - 1) && absY == (ny - 1) && absZ == (nz - 1)) {

                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == 0 && absY == 0 && absZ == 0) {

                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);
                                }

                                else if (absX == 0 && absY == (ny - 1) && absZ == (nz - 1)) {

                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);

                                }
                                else if (absX == (nx - 1) && absY == 0 && absZ == 0) {

                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(1);

                                }
                                else if (absX == (nx - 1) && absY == ny - 1 && absZ == nz - 1) {

                                    delXYZ.push_back(-1); delXYZ.push_back(0); delXYZ.push_back(0); delXYZ.push_back(-1);
                                }
                                else {

                                    delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1); delXYZ.push_back(0);
                                    delXYZ.push_back(0); delXYZ.push_back(1); delXYZ.push_back(0); delXYZ.push_back(-1);
                                }
                                std::vector<plint> del_iX, del_iY, del_iZ;
                                for (plint iT = 0; iT < nbrs; ++iT) {
                                    plint delx = delXYZ[iT * 2], dely = delXYZ[iT * 2 + 1], delz = delXYZ[iT * 2 + 1];
                                    if (util::roundToInt(lattices[distLloc]->get(iXd + delx, iYd + dely, iZd + delz).computeDensity()) > 0) {
                                        del_iX.push_back(delx); del_iY.push_back(dely); del_iZ.push_back(delz);
                                    }
                                }
                                T newAge = 0; bool updateflag = 1;
                                plint age = util::roundToInt(lattices[ageLloc]->get(iXa, iYa, iZa).computeDensity());
                                if (age == 0) { newAge = 1; }
                                else if (age == 1 && (B - Bmax) > -thrd) {
                                    for (size_t iT = 0; iT < del_iX.size(); ++iT) {
                                        plint delx = del_iX[iT], dely = del_iY[iT], delz = del_iZ[iT];
                                        plint nbrAge = util::roundToInt(lattices[ageLloc]->get(iXa + delx, iYa + dely, iZa + delz).computeDensity());
                                        plint nbrDist = util::roundToInt(lattices[distLloc]->get(iXd + delx, iYd + dely, iZd + delz).computeDensity());
                                        if (nbrDist > 0 && nbrAge == 0) { updateflag = 0; break; }
                                    }
                                    if (updateflag == 1) { newAge = 2; }
                                }
                                else {
                                    plint count = 0;
                                    for (size_t iT = 0; iT < del_iX.size(); ++iT) {
                                        plint delx = del_iX[iT], dely = del_iY[iT], delz = del_iZ[iT];
                                        plint nbrAge = util::roundToInt(lattices[ageLloc]->get(iXa + delx, iYa + dely, iZa + delz).computeDensity());
                                        T nbrB = lattices[bMtLloc]->get(iXb + delx, iYb + dely, iZb + delz).computeDensity();
                                        if (nbrAge >= age && (nbrB - Bmax) > -thrd) { ++count; }
                                    }
                                    if (count == del_iX.size()) { newAge = (T)(age + 1); }
                                    else { updateflag = 0; }
                                }
                                if (updateflag == 1) {
                                    Array<T, 7> g;
                                    g[0] = (T)(newAge - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(newAge - 1) / 8;
                                    lattices[ageLloc]->get(iXa, iYa, iZa).setPopulations(g);
                                }
                            }
                        }
                    }
                }
            }
                virtual BlockDomain::DomainT appliesTo() const {
                    return BlockDomain::bulk;
                }
                virtual updateAgeDistance3D<T, Descriptor>* clone() const {
                    return new updateAgeDistance3D<T, Descriptor>(*this);
                }
                void getTypeOfModification(std::vector<modif::ModifT>&modified) const {
                    modified[0] = modif::staticVariables;
                    modified[1] = modif::nothing;
                    modified[2] = modif::nothing;
                }
        private:
            T Bmax;
            plint nx, ny, nz;
        };



            // create a domain distance scalarfield3d
        
            template<typename T1>
            class createDistanceDomain : public BoxProcessingFunctional3D_S<T1>
            {
            public:
                createDistanceDomain(std::vector<std::vector<std::vector<plint>>> distVec_) : distVec(distVec_)
                {}
                virtual void process(Box3D domain, ScalarField3D<T1>& field) {
                
                    Dot3D absoluteOffset = field.getLocation();
                  
                    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                        plint absX = iX + absoluteOffset.x;
                        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                           plint absY = iY + absoluteOffset.y;                           
                           for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                                plint absZ = iZ + absoluteOffset.z;
                                field.get(iX,iY,iZ) = distVec[absX][absY][absZ];
                           } 
                        }
                    }
                }
                virtual createDistanceDomain<T1>* clone() const {
                    return new createDistanceDomain<T1>(*this);
                }
                virtual BlockDomain::DomainT appliesTo() const {
                    return BlockDomain::bulkAndEnvelope;
                }
                void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                    modified[0] = modif::staticVariables;
                }
            private:
                std::vector<std::vector<std::vector<plint>>> distVec;
            }; 

            // create a domain distance scalarfield3d
            template<typename T1>
            class createAgeDomain : public BoxProcessingFunctional3D_S<T1>
            {
            public:
                createAgeDomain(std::vector<plint> pore_, plint bb_, plint solid_) : pore(pore_), bb(bb_), solid(solid_)
                {}
                virtual void process(Box3D domain, ScalarField3D<T1>& field) {
                    for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
                        for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                            for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                                plint mask = field.get(iX, iY, iZ);
                                if (mask == solid || mask == bb) { field.get(iX, iY, iZ) = -1; }
                                else {
                                    bool poreflag = 0;
                                    for (size_t iP = 0; iP < pore.size(); ++iP) { if (mask == pore[iP]) { poreflag = 1; break; } }
                                    if (poreflag == 1) { field.get(iX, iY, iZ) = 0; }
                                    else { field.get(iX, iY, iZ) = 1; }
                                }
                            }
                        }
                    }
                }
                virtual createAgeDomain<T1>* clone() const {
                    return new createAgeDomain<T1>(*this);
                }
                virtual BlockDomain::DomainT appliesTo() const {
                    return BlockDomain::bulkAndEnvelope;
                }
                void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
                    modified[0] = modif::staticVariables;
                }
            private:
                std::vector<plint> pore;
                plint bb, solid;
            };

            /* ===============================================================================================================
               ========================================== REDUCTIVE DATA PROCESSORS ==========================================
               =============================================================================================================== */




// Data reduction is a process that reduced the volume of original data and represents it in a much smaller volume. 

template<typename T1>
class MaskedBoxScalarCountFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T1>
{
public:
    MaskedBoxScalarCountFunctional3D(plint mask_) : countId(this->getStatistics().subscribeSum()), mask(mask_)
    {}
    virtual void process(Box3D domain, ScalarField3D<T1>& scalar) {
        BlockStatistics& statistics = this->getStatistics();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint tmpMask = util::roundToInt(scalar.get(iX, iY, iZ));
                    if (tmpMask == mask) {
                        statistics.gatherSum(countId, (int)1);

                    }

                }
            }
        }
    }
    virtual MaskedBoxScalarCountFunctional3D<T1>* clone() const {
        return new MaskedBoxScalarCountFunctional3D<T1>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    plint getCount() const {
        // The sum is internally computed on floating-point values.
        // If T is integer, the value must be rounded at the end.
        plint doubleSum = this->getStatistics().getSum(countId);
        // if (std::numeric_limits<T>::is_integer) {
        //     return (T) util::roundToInt(doubleSum);
        // }
        return (T)doubleSum;
    }
private:
    plint countId;
    plint mask;
};

template<typename T1>
plint MaskedScalarCounts(Box3D domain, MultiScalarField3D<T1>& field, plint mask) {
    MaskedBoxScalarCountFunctional3D<T1> functional = MaskedBoxScalarCountFunctional3D<T1>(mask);
    applyProcessingFunctional(functional, domain, field);
    return functional.getCount();
}

// calculate RMSE for convergence checking
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class BoxLatticeRMSEFunctional3D : public ReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>
{
public:
    BoxLatticeRMSEFunctional3D() : sumId(this->getStatistics().subscribeSum())
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice0, BlockLattice3D<T2, Descriptor2>& lattice1) {
        BlockStatistics& statistics = this->getStatistics();
        Dot3D offset_01 = computeRelativeDisplacement(lattice0, lattice1);
        for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
            for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                    plint iX1 = iX0 + offset_01.x; plint iY1 = iY0 + offset_01.y; plint iZ1 = iZ0 + offset_01.z;
                    T deltaC = lattice0.get(iX0, iY0, iZ0).computeDensity() - lattice1.get(iX1, iY1, iZ1).computeDensity();
                    T RMSE = deltaC * deltaC;
                    statistics.gatherSum(sumId, RMSE);
                }
            }
        }
    }
    virtual BoxLatticeRMSEFunctional3D<T1, Descriptor1, T2, Descriptor2>* clone() const {
        return new BoxLatticeRMSEFunctional3D<T1, Descriptor1, T2, Descriptor2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getCount() const {
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        double doubleSum = this->getStatistics().getSum(sumId);
        if (std::numeric_limits<T>::is_integer) {
            return (T)util::roundToInt(doubleSum);
        }
        return (T)doubleSum;
    }
private:
    plint sumId;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
T computeRMSE(Box3D domain, MultiBlockLattice3D<T1, Descriptor1>& lattice0, MultiBlockLattice3D<T2, Descriptor2>& lattice1, T poreLen) {
    BoxLatticeRMSEFunctional3D<T1, Descriptor1, T2, Descriptor2> functional;
    applyProcessingFunctional(functional, domain, lattice0, lattice1);
    return std::sqrt(functional.getCount() / poreLen);
}
