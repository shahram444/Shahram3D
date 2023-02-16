// nonlocal biomass reaction
template<typename T, template<typename U> class Descriptor>
class localReactionLattices3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    localReactionLattices3D(T dt_) : dt(dt_)
    {}
    // lattices[0] = solute 1 concentration
    // lattices[1] = surface-associated biomass 1
    // lattices[2] = mask field lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        Dot3D offset_01 = computeRelativeDisplacement(*lattices[0], *lattices[1]);
        Dot3D offset_02 = computeRelativeDisplacement(*lattices[0], *lattices[2]);
        for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
            plint iX2 = iX0 + offset_02.x;
            for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                plint iY2 = iY0 + offset_02.y;
                for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                    plint iZ2 = iZ0 + offset_02.z;
                    plint mask = util::roundToInt(lattices[2]->get(iX2, iY2, iZ2).computeDensity());
                    if ((mask > 1 && mask < 5) || mask > 6) { // biomass
                        Array<T, 7> g;
                        plint iX1 = iX0 + offset_01.x; plint iY1 = iY0 + offset_01.y; plint iZ1 = iZ0 + offset_01.z;

                        T c0 = lattices[0]->get(iX0, iY0, iZ0).computeDensity(); // (mM)
                        T B2 = lattices[1]->get(iX1, iY1, iZ0).computeDensity(); // sessile (gdw/L)

                        call_biomassGrowth(c0, B2, dt);

                        lattices[0]->get(iX0, iY0, iZ0).getPopulations(g);
                        g[0] += (T)c0 / 4; g[1] += (T)c0 / 8; g[2] += (T)c0 / 8; g[3] += (T)c0 / 8; g[4] += (T)c0 / 8; g[5] += (T)c0 / 8; g[6] += (T)c0 / 8;
                        lattices[0]->get(iX0, iY0, iZ0).setPopulations(g);

                        lattices[1]->get(iX1, iY1, iZ1).getPopulations(g);
                        g[0] += (T)B2 / 4; g[1] += (T)B2 / 8; g[2] += (T)B2 / 8; g[3] += (T)B2 / 8; g[4] += (T)B2 / 8; g[5] += (T)c0 / 8; g[6] += (T)c0 / 8;
                        lattices[1]->get(iX1, iY1, iZ1).setPopulations(g);
                    }
                }
            }
        }

    }
   
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual localReactionLattices3D<T, Descriptor>* clone() const {
        return new localReactionLattices3D<T, Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::staticVariables;
        modified[2] = modif::nothing;
    }
private:
    T dt;
};

// nonlocal biomass reaction
template<typename T, template<typename U> class Descriptor>
class updateTOTbMassLattices3D : public LatticeBoxProcessingFunctional3D<T, Descriptor>
{
public:
    updateTOTbMassLattices3D(plint length_) : length(length_)
    {}
    // lattices[0~(#ofbM-1)] = original biomass lattices
    // lattices[#ofbM~(len-3)] = copy biomass lattices
    // lattices[len-2] = total biomass lattice
    // lattices[len-1] = mask lattice
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T, Descriptor>*> lattices) {
        std::vector<Dot3D> vec_offset;
        for (plint iL = 0; iL < length; ++iL) {
            vec_offset.push_back(computeRelativeDisplacement(*lattices[0], *lattices[iL]));
        }
        plint numbM = (length - 2) / 2;
        for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
            for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                    plint iXm = iX0 + vec_offset[length - 1].x; plint iYm = iY0 + vec_offset[length - 1].y; plint iZm = iZ0 + vec_offset[length - 1].z;
                    plint mask = util::roundToInt(lattices[length - 1]->get(iXm, iYm, iZm).computeDensity());
                    if (mask > 2) {
                        T bmass = 0;
                        for (plint iB = 0; iB < numbM; ++iB) {
                            plint iXb = iX0 + vec_offset[iB].x; plint iYb = iY0 + vec_offset[iB].y; plint iZb = iZ0 + vec_offset[iB].z;
                            bmass += lattices[iB]->get(iXb, iYb, iZb).computeDensity();
                        }
                        Array<T, 7> g;
                        g[0] = (T)(bmass - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(bmass - 1) / 8;
                        plint iXt = iX0 + vec_offset[length - 2].x; plint iYt = iY0 + vec_offset[length - 2].y; plint iZt = iZ0 + vec_offset[length - 2].z;
                        lattices[length - 2]->get(iXt, iYt, iZt).setPopulations(g);
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    virtual updateTOTbMassLattices3D<T, Descriptor>* clone() const {
        return new updateTOTbMassLattices3D<T, Descriptor>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (plint iL = 0; iL < length; ++iL) {
            if (iL == length - 2) {
                modified[iL] = modif::staticVariables;
            }
            else {
                modified[iL] = modif::nothing;
            }
        }
    }
private:
    plint length;
};

// Update local dynamics of the mask lattice
template<typename T1, template<typename U> class Descriptor1, typename T2, template<typename U> class Descriptor2>
class updateMaskLatticeDynamics3D : public BoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>
{
public:
    updateMaskLatticeDynamics3D(T bMf_, T Bmax_, plint pore_, plint solid_, plint bb_)
        : bMthrd(bMf_* Bmax_), pore(pore_), solid(solid_), bb(bb_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice1, BlockLattice3D<T2, Descriptor2>& lattice2) {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint mask = util::roundToInt(lattice1.get(iX, iY, iZ).computeDensity());
                    if (mask != pore && mask != solid && mask != bb) {
                        T bmass = lattice2.get(iX, iY, iZ).computeDensity();
                        plint omega = util::roundToInt(lattice2.get(iX, iY, iZ).getDynamics().getOmega());
                        if (bmass >= bMthrd && omega == 0) {
                            lattice2.get(iX, iY, iZ).getDynamics().setOmega(1.);
                        }

                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    virtual updateMaskLatticeDynamics3D<T1, Descriptor1, T2, Descriptor2>* clone() const {
        return new updateMaskLatticeDynamics3D<T1, Descriptor1, T2, Descriptor2>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;
    }
private:
    T bMthrd;
    plint pore, solid, bb;
};


// lattice divisions for MPI


template<typename T1, template<typename U> class Descriptor, typename T2>
class latticeXYZ3D : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
{
public:
    latticeXYZ3D(bool biId_) : biId(biId_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint iX1 = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint iY1 = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint iZ1 = iZ + offset.z;
                    if (biId == 0) {
                        field.get(iX1, iY1, iZ1) = iX;
                    }
                    else if (biId == 1) {
                        field.get(iX1, iY1, iZ1) = iY;

                    }

                    else {

                        field.get(iX1, iY1, iZ1) = iZ;

                    }
                  
                }
            }
        }
    }
    virtual latticeXYZ3D<T1, Descriptor, T2>* clone() const {
        return new latticeXYZ3D<T1, Descriptor, T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    bool biId;
};

//seconde version of lattice divisions for MPI
// template<typename T1, template<typename U> class Descriptor, typename T2>
// class latticeXY3D : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
// {
// public:
//     latticeXY3D(bool biId_) : biId(biId_)
//     {}
//     virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
//         Dot3D offset = computeRelativeDisplacement(lattice, field);
//         for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
//             plint iX1 = iX + offset.x;
//             for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
//                 plint iY1 = iY + offset.y;
//                 for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
//                     plint iZ1 = iZ + offset.z;
//                     if (biId == 0) {
//                         field.get(iX1, iY1, iZ1) = iX;
//                     }
//                     else if (biId == 1) {
//                         field.get(iX1, iY1, iZ1) = iY;
//                     }
//                     else {
//                         field.get(iX1, iY1, iZ1) = iZ;
//                     }
//                 }
//             }
//         }
//     }
//     virtual latticeXY3D<T1, Descriptor, T2>* clone() const {
//         return new latticeXY3D<T1, Descriptor, T2>(*this);
//     }
//     virtual BlockDomain::DomainT appliesTo() const {
//         return BlockDomain::bulkAndEnvelope;
//     }
//     void getTypeOfModification(std::vectormodif::ModifT& modified) const {
//         modified[0] = modif::nothing;
//         modified[1] = modif::staticVariables;
//     }
// private:
//     bool biId;
// };


// absolute lattice indices
template<typename T1, template<typename U> class Descriptor, typename T2>
class latticeAbsoluteXYZ3D : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
{
public:
    latticeAbsoluteXYZ3D(bool biId_) : biId(biId_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor>& lattice, ScalarField3D<T2>& field) {
        Dot3D offset = computeRelativeDisplacement(lattice, field);
        Dot3D absoluteOffset = lattice.getLocation();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            plint absX = iX + absoluteOffset.x;
            plint iX1 = iX + offset.x;
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                plint absY = iY + absoluteOffset.y;
                plint iY1 = iY + offset.y;
                for (plint iZ = domain.z0; iZ <= domain.y1; ++iZ) {
                    plint absZ = iZ + absoluteOffset.z;
                    plint iZ1 = iZ + offset.z;
                    if (biId == 0) {
                        field.get(iX1, iY1,iZ1) = absX;
                    }
                    
                    else if (biId == 1) {
                    field.get(iX1, iY1, iZ1) = absY;
                   
                    }

                          else { 
                          field.get(iX1, iY1, iZ1) = absZ;
                          } 
                }
            }
        }
    }
    virtual latticeAbsoluteXYZ3D<T1, Descriptor, T2>* clone() const {
        return new latticeAbsoluteXYZ3D<T1, Descriptor, T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    bool biId;
};

template<typename T1, template<typename U> class Descriptor, typename T2>
class initializeMaskedScalarLattice : public BoxProcessingFunctional3D_LS<T1, Descriptor, T2>
{
public:
    initializeMaskedScalarLattice(T rho0_, plint id_) : rho0(rho0_), id(id_)
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
                    if (mask == id) {
                        Array<T, 7> g;
                        g[0] = (T)(rho0 - 1) / 4; g[1] = g[2] = g[3] = g[4] = g[5] = g[6] = (T)(rho0 - 1) / 8;
                        lattice.get(iX, iY, iZ).setPopulations(g);
                    }
                }
            }
        }
    }
    virtual initializeMaskedScalarLattice<T1, Descriptor, T2>* clone() const {
        return new initializeMaskedScalarLattice<T1, Descriptor, T2>(*this);
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
private:
    T rho0;
    plint id;
};


/* ===============================================================================================================
   ========================================== REDUCTIVE DATA PROCESSORS ==========================================
   =============================================================================================================== */

   // count the number of a certain mask
template<typename T1, template<typename U1> class Descriptor1>
class BoxScalarSelectedSumFunctional3D : public ReductiveBoxProcessingFunctional3D_L<T1, Descriptor1>
{
public:
    BoxScalarSelectedSumFunctional3D(plint mask_) : countId(this->getStatistics().subscribeSum()), mask(mask_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice) {
        BlockStatistics& statistics = this->getStatistics();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint tmpMask = util::roundToInt(lattice.get(iX, iY, iZ).computeDensity());
                    if (tmpMask == mask) {
                        statistics.gatherSum(countId, (int)1);
                    }
                }
            }
        }
    }
    virtual BoxScalarSelectedSumFunctional3D<T1, Descriptor1>* clone() const {
        return new BoxScalarSelectedSumFunctional3D<T1, Descriptor1>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getCount() const {
        // The sum is internally computed on floating-point values.
        // If T is integer, the value must be rounded at the end.
        double doubleSum = this->getStatistics().getSum(countId);
        if (std::numeric_limits<T>::is_integer) {
            return (T)util::roundToInt(doubleSum);
        }
        return (T)doubleSum;
    }
private:
    plint countId;
    plint mask;
};

template<typename T1, template<typename U1> class Descriptor1>
T countLatticeMaskNumbers(Box3D domain, MultiBlockLattice3D<T1, Descriptor1>& lattice, plint mask) {
    BoxScalarSelectedSumFunctional3D<T1, Descriptor1> functional = BoxScalarSelectedSumFunctional3D<T1, Descriptor1>(mask);
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getCount();
}

// sum up lattice omegas
template<typename T1, template<typename U1> class Descriptor1>
class SumLatticeCellOmegas3D : public ReductiveBoxProcessingFunctional3D_L<T1, Descriptor1>
{
public:
    SumLatticeCellOmegas3D() : omegasum(this->getStatistics().subscribeSum())
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice) {
        BlockStatistics& statistics = this->getStatistics();
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY) {
                for (plint iZ = domain.z0; iZ <= domain.z1; ++iZ) {
                    plint omega = util::roundToInt(lattice.get(iX, iY, iZ).getDynamics().getOmega());
                    statistics.gatherSum(omegasum, omega);
                }
            }
        }
    }
    virtual SumLatticeCellOmegas3D<T1, Descriptor1>* clone() const {
        return new SumLatticeCellOmegas3D<T1, Descriptor1>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getSum() const {
        // The sum is internally computed on floating-point values.
        // If T is integer, the value must be rounded at the end.
        double doubleSum = this->getStatistics().getSum(omegasum);
        if (std::numeric_limits<T>::is_integer) {
            return (T)util::roundToInt(doubleSum);
        }
        return (T)doubleSum;
    }
private:
    plint omegasum;
};

template<typename T1, template<typename U1> class Descriptor1>
T sumUpLatticeOmegas(Box3D domain, MultiBlockLattice3D<T1, Descriptor1>& lattice) {
    SumLatticeCellOmegas3D<T1, Descriptor1> functional = SumLatticeCellOmegas3D<T1, Descriptor1>();
    applyProcessingFunctional(functional, domain, lattice);
    return functional.getSum();
}

// calculate the maximum density with a mask lattice
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class MaskedBoxLatticeMaxFunctional3D : public ReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>
{
public:
    MaskedBoxLatticeMaxFunctional3D(plint pore_, plint solid_, plint bb_) : maxLatticeId(this->getStatistics().subscribeMax()), pore(pore_), solid(solid_), bb(bb_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice0, BlockLattice3D<T2, Descriptor2>& lattice1) {
        BlockStatistics& statistics = this->getStatistics();
        Dot3D offset_01 = computeRelativeDisplacement(lattice0, lattice1);
        for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
            for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                    plint mask = util::roundToInt(lattice1.get(iX0 + offset_01.x, iY0 + offset_01.y, iZ0 + offset_01.z).computeDensity());
                    // if (mask > 2) {
                    if (mask != pore && mask != solid && mask != bb) {
                        T max = lattice0.get(iX0, iY0, iZ0).computeDensity();
                        statistics.gatherMax(maxLatticeId, max);
                    }
                }
            }
        }
    }
    virtual MaskedBoxLatticeMaxFunctional3D<T1, Descriptor1, T2, Descriptor2>* clone() const {
        return new MaskedBoxLatticeMaxFunctional3D<T1, Descriptor1, T2, Descriptor2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMaxLattice() const {
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        double doubleMax = this->getStatistics().getMax(maxLatticeId);
        if (std::numeric_limits<T>::is_integer) {
            return (T)util::roundToInt(doubleMax);
        }
        return (T)doubleMax;
    }
private:
    plint maxLatticeId, pore, solid, bb;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
T computeMaskedLatticeMax(Box3D domain, MultiBlockLattice3D<T1, Descriptor1>& lattice0, MultiBlockLattice3D<T2, Descriptor2>& lattice1, plint pore, plint solid, plint bb) {
    MaskedBoxLatticeMaxFunctional3D<T1, Descriptor1, T2, Descriptor2> functional = MaskedBoxLatticeMaxFunctional3D<T1, Descriptor1, T2, Descriptor2>(pore, solid, bb);
    applyProcessingFunctional(functional, domain, lattice0, lattice1);
    return functional.getMaxLattice();
}

// calculate the minimum density with a mask lattice
template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
class MaskedBoxLatticeMinFunctional3D : public ReductiveBoxProcessingFunctional3D_LL<T1, Descriptor1, T2, Descriptor2>
{
public:
    MaskedBoxLatticeMinFunctional3D(plint mask1_, plint mask2_) : maxLatticeId(this->getStatistics().subscribeMax()), mask1(mask1_), mask2(mask2_)
    {}
    virtual void process(Box3D domain, BlockLattice3D<T1, Descriptor1>& lattice0, BlockLattice3D<T2, Descriptor2>& lattice1) {
        BlockStatistics& statistics = this->getStatistics();
        Dot3D offset_01 = computeRelativeDisplacement(lattice0, lattice1);
        for (plint iX0 = domain.x0; iX0 <= domain.x1; ++iX0) {
            for (plint iY0 = domain.y0; iY0 <= domain.y1; ++iY0) {
                for (plint iZ0 = domain.z0; iZ0 <= domain.z1; ++iZ0) {
                    plint tmpMask = lattice1.get(iX0 + offset_01.x, iY0 + offset_01.y, iZ0 + offset_01.z).computeDensity();
                    if (tmpMask == mask1 || tmpMask == mask2) {
                        T min = -lattice0.get(iX0, iY0, iZ0).computeDensity();
                        statistics.gatherMax(maxLatticeId, min);
                    }
                }
            }
        }
    }
    virtual MaskedBoxLatticeMinFunctional3D<T1, Descriptor1, T2, Descriptor2>* clone() const {
        return new MaskedBoxLatticeMinFunctional3D<T1, Descriptor1, T2, Descriptor2>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getMinLattice() const {
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        double doubleMin = -this->getStatistics().getMax(maxLatticeId);
        if (std::numeric_limits<T>::is_integer) {
            return (T)util::roundToInt(doubleMin);
        }
        return (T)doubleMin;
    }
private:
    plint maxLatticeId;
    plint mask1;
    plint mask2;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2>
T computeMaskedLatticeMin(Box3D domain, MultiBlockLattice3D<T1, Descriptor1>& lattice0, MultiBlockLattice3D<T2, Descriptor2>& lattice1, plint mask1, plint mask2) {
    MaskedBoxLatticeMinFunctional3D<T1, Descriptor1, T2, Descriptor2> functional = MaskedBoxLatticeMinFunctional3D<T1, Descriptor1, T2, Descriptor2>(mask1, mask2);
    applyProcessingFunctional(functional, domain, lattice0, lattice1);
    return functional.getMinLattice();
}

