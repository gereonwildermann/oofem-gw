/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "sm/Elements/Bars/rebar3d.h"
#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "engngm.h"
#include "mathfem.h"
#include "parametermanager.h"
#include "paramkey.h"
#include "fei3dlinelin.h"
#include "node.h"
#include "gaussintegrationrule.h"
#include "classfactory.h"
#include <iostream>

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
ParamKey Rebar3d::IPK_Rebar3d_nlgeoflag("nlgeo");

REGISTER_Element(Rebar3d);

FEI3dLineLin Rebar3d::interp;

Rebar3d::Rebar3d(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    nlGeometry = 0; // Geometrical nonlinearities disabled as default
    numberOfDofMans = 2;
}


FEInterpolation *Rebar3d::giveInterpolation() const { return & interp; }


Interface *
Rebar3d::giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    }

    //OOFEM_LOG_INFO("Interface on Rebar3d element not supported");
    return NULL;
}

double
Rebar3d::computeCurrentVolume(TimeStep *tStep)
{
    double vol = 0.0;
    for ( auto &gp : * this->giveDefaultIntegrationRulePtr() ) {
        FloatArray F;
        FloatMatrix Fm;


        computeDeformationGradientVector(F, gp, tStep);
        Fm.beMatrixForm(F);
        double J = Fm.giveDeterminant();

        FEInterpolation *interpolation = this->giveInterpolation();
        double detJ = fabs( ( interpolation->giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) ) );

        vol += gp->giveWeight() * detJ * J;
    }

    return vol;
}


void
Rebar3d::computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf({ D_u, D_v, D_w }, VM_Total, tStep, u); // solution vector
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt(gp, B);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}


void
Rebar3d::computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep);
    answer = cs->giveFirstPKStresses(vF, gp, tStep);
}

void
Rebar3d::computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = this->giveStructuralCrossSection();

    FloatArray vF;
    this->computeDeformationGradientVector(vF, gp, tStep);
    cs->giveCauchyStresses(answer, gp, vF, tStep);
}

void
Rebar3d::giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix B;
    FloatArray vStress, vStrain, u;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( auto &gp : * this->giveDefaultIntegrationRulePtr() ) {
        StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( this->giveCrossSection()->giveMaterial(gp)->giveStatus(gp) );

        // Engineering (small strain) stress
        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B);
            if ( useUpdatedGpRecord == 1 ) {
                vStress = matStat->giveStressVector();
            } else {
                ///@todo Is this really what we should do for inactive elements?
                if ( !this->isActivated(tStep) ) {
                    vStrain.resize(StructuralMaterial::giveSizeOfVoigtSymVector(gp->giveMaterialMode() ) );
                    vStrain.zero();
                }
                vStrain.beProductOf(B, u);
                this->computeStressVector(vStress, vStrain, gp, tStep);
            }
        } else if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->giveCVector();
                } else {
                    this->computeCauchyStressVector(vStress, gp, tStep);
                }

                this->computeBmatrixAt(gp, B);
            } else { // First Piola-Kirchhoff stress
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->givePVector();
                } else {
                    this->computeFirstPKStressVector(vStress, gp, tStep);
                    ///@todo This is actaully inefficient since it constructs B and twice and collects the nodal unknowns over and over.
                }

                this->computeBHmatrixAt(gp, B);
            }
        }

        if ( vStress.giveSize() == 0 ) { /// @todo is this check really necessary?
            break;
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = this->computeVolumeAround(gp, tStep);

        if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            if ( vStress.giveSize() == 9 ) {
                FloatArray stressTemp;
                StructuralMaterial::giveReducedVectorForm(stressTemp, vStress, gp->giveMaterialMode() );
                answer.plusProduct(B, stressTemp, dV);
            } else {
                answer.plusProduct(B, vStress, dV);
            }
        } else {
            if ( vStress.giveSize() == 6 ) {
                // It may happen that e.g. plane strain is computed
                // using the default 3D implementation. If so,
                // the stress needs to be reduced.
                // (Note that no reduction will take place if
                //  the simulation is actually 3D.)
                FloatArray stressTemp;
                StructuralMaterial::giveReducedSymVectorForm(stressTemp, vStress, gp->giveMaterialMode() );
                answer.plusProduct(B, stressTemp, dV);
            } else {
                answer.plusProduct(B, vStress, dV);
            }
        }
    }

    // If inactive: update fields but do not give any contribution to the internal forces
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}


void
Rebar3d::giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer,
                                                                   TimeStep *tStep, int useUpdatedGpRecord)
{
    /*
     * Returns nodal representation of real internal forces computed from first Piola-Kirchoff stress
     * if useGpRecord == 1 then stresses stored in the gp are used, otherwise stresses are computed
     * this must be done if you want internal forces after element->updateYourself() has been called
     * for the same time step.
     * The integration procedure uses an integrationRulesArray for numerical integration.
     * Each integration rule is considered to represent a separate sub-cell/element. Typically this would be used when
     * integration of the element domain needs special treatment, e.g. when using the XFEM.
     */

    FloatMatrix B;
    FloatArray vStress, vStrain;

    IntArray irlocnum;
    FloatArray *m = & answer, temp;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();


    // loop over individual integration rules
    for ( auto &iRule : integrationRulesArray ) {
        for ( GaussPoint *gp : * iRule ) {
            StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );

            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                if ( useUpdatedGpRecord == 1 ) {
                    vStress = matStat->giveStressVector();
                } else {
                    this->computeStrainVector(vStrain, gp, tStep);
                    this->computeStressVector(vStress, vStrain, gp, tStep);
                }
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Cauchy stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = matStat->giveCVector();
                    } else {
                        this->computeCauchyStressVector(vStress, gp, tStep);
                    }

                    this->computeBmatrixAt(gp, B);
                } else { // First Piola-Kirchhoff stress
                    if ( useUpdatedGpRecord == 1 ) {
                        vStress = matStat->givePVector();
                    } else {
                        this->computeFirstPKStressVector(vStress, gp, tStep);
                    }

                    this->computeBHmatrixAt(gp, B);
                }
            }

            if ( vStress.giveSize() == 0 ) { //@todo is this really necessary?
                break;
            }

            // compute nodal representation of internal forces at nodes as f = B^T*stress dV
            double dV = this->computeVolumeAround(gp, tStep);
            m->plusProduct(B, vStress, dV);

            // localize irule contribution into element matrix
            if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, * iRule) ) {
                answer.assemble(* m, irlocnum);
                m->clear();
            }
        }
    }

    // if inactive: update fields but do not give any contribution to the structure
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}






void
Rebar3d::computeStiffnessMatrix(FloatMatrix &answer,
                                            MatResponseMode rMode, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();

    if ( !this->isActivated(tStep) ) {
        return;
    }

    // Compute matrix from material stiffness (total stiffness for small def.) - B^T * dS/dE * B
    if ( integrationRulesArray.size() == 1 ) {
        FloatMatrix B, D, DB;
        for ( auto &gp : * this->giveDefaultIntegrationRulePtr() ) {
            // Engineering (small strain) stiffness
            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B);
                    /// @todo We probably need overloaded function (like above) here as well.
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B);
                    this->computeConstitutiveMatrix_dPdF_At(D, rMode, gp, tStep);
                }
            }

            double dV = this->computeVolumeAround(gp, tStep);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                answer.plusProductSymmUpper(B, DB, dV);
            } else {
                answer.plusProductUnsym(B, DB, dV);
            }
        }

        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
            FloatMatrix initialStressMatrix;
            this->computeInitialStressMatrix(initialStressMatrix, tStep);
            answer.add(initialStressMatrix);
        }
    } else { /// @todo Verify that it works with large deformations
        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
            OOFEM_ERROR("Updated lagrangian not supported yet");
        }

        int iStartIndx, iEndIndx, jStartIndx, jEndIndx;
        FloatMatrix Bi, Bj, D, Dij, DBj;
        for ( int i = 0; i < ( int ) integrationRulesArray.size(); i++ ) {
            iStartIndx = integrationRulesArray [ i ]->getStartIndexOfLocalStrainWhereApply();
            iEndIndx   = integrationRulesArray [ i ]->getEndIndexOfLocalStrainWhereApply();
            for ( int j = 0; j < ( int ) integrationRulesArray.size(); j++ ) {
                IntegrationRule *iRule;
                jStartIndx = integrationRulesArray [ j ]->getStartIndexOfLocalStrainWhereApply();
                jEndIndx   = integrationRulesArray [ j ]->getEndIndexOfLocalStrainWhereApply();
                if ( i == j ) {
                    iRule = integrationRulesArray [ i ].get();
                } else if ( integrationRulesArray [ i ]->giveNumberOfIntegrationPoints() < integrationRulesArray [ j ]->giveNumberOfIntegrationPoints() ) {
                    iRule = integrationRulesArray [ i ].get();
                } else {
                    iRule = integrationRulesArray [ j ].get();
                }

                for ( GaussPoint *gp : * iRule ) {
                    // Engineering (small strain) stiffness dSig/dEps
                    if ( nlGeometry == 0 ) {
                        this->computeBmatrixAt(gp, Bi, iStartIndx, iEndIndx);
                        this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
                    } else if ( nlGeometry == 1 ) {
                        if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                            this->computeBmatrixAt(gp, Bi);
                            cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                        } else { // Material stiffness dP/dF
                            this->computeBHmatrixAt(gp, Bi);
                            this->computeConstitutiveMatrix_dPdF_At(D, rMode, gp, tStep);
                            //cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                        }
                    }


                    if ( i != j ) {
                        if ( nlGeometry == 0 ) {
                            this->computeBmatrixAt(gp, Bj, jStartIndx, jEndIndx);
                        } else if ( nlGeometry == 1 ) {
                            if ( this->domain->giveEngngModel()->giveFormulation() == AL ) {
                                this->computeBmatrixAt(gp, Bj);
                            } else {
                                this->computeBHmatrixAt(gp, Bj);
                            }
                        }
                    } else {
                        Bj  = Bi;
                    }

                    Dij.beSubMatrixOf(D, iStartIndx, iEndIndx, jStartIndx, jEndIndx);
                    double dV = this->computeVolumeAround(gp, tStep);
                    DBj.beProductOf(Dij, Bj);
                    if ( matStiffSymmFlag ) {
                        answer.plusProductSymmUpper(Bi, DBj, dV);
                    } else {
                        answer.plusProductUnsym(Bi, DBj, dV);
                    }
                }
            }
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
Rebar3d::computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer,
                                                                 MatResponseMode rMode, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    bool matStiffSymmFlag = cs->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    FloatMatrix temp;
    FloatMatrix *m = & answer;
    if ( this->giveInterpolation() && this->giveInterpolation()->hasSubPatchFormulation() ) {
        m = & temp;
    }

    // Compute matrix from material stiffness
    FloatMatrix B, D, DB;
    IntArray irlocnum;
    for ( auto &iRule : integrationRulesArray ) {
        for ( auto &gp : * iRule ) {
            if ( nlGeometry == 0 ) {
                this->computeBmatrixAt(gp, B);
                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
            } else if ( nlGeometry == 1 ) {
                if ( this->domain->giveEngngModel()->giveFormulation() == AL ) { // Material stiffness dC/de
                    this->computeBmatrixAt(gp, B);
                    cs->giveStiffnessMatrix_dCde(D, rMode, gp, tStep);
                } else { // Material stiffness dP/dF
                    this->computeBHmatrixAt(gp, B);
                    this->computeConstitutiveMatrix_dPdF_At(D, rMode, gp, tStep);
                    //cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);
                }
            }

            double dV = this->computeVolumeAround(gp, tStep);
            DB.beProductOf(D, B);
            if ( matStiffSymmFlag ) {
                m->plusProductSymmUpper(B, DB, dV);
            } else {
                m->plusProductUnsym(B, DB, dV);
            }
        }

        // localize irule contribution into element matrix
        if ( this->giveIntegrationRuleLocalCodeNumbers(irlocnum, * iRule) ) {
            answer.assemble(* m, irlocnum);
            m->clear();
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}

void
Rebar3d::initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm = this->domain->elementPPM;
    StructuralElement::initializeFrom(ir, priority);
    PM_UPDATE_PARAMETER(nlGeometry, ppm, ir, this->number, IPK_Rebar3d_nlgeoflag, priority) ;
}

void Rebar3d::giveInputRecord(DynamicInputRecord &input)
{
    StructuralElement::giveInputRecord(input);

    input.setField(nlGeometry, IPK_Rebar3d_nlgeoflag.getNameCStr());
}

int
Rebar3d::checkConsistency()
{
    if ( this->nlGeometry == 2 ) {
        OOFEM_ERROR("nlGeometry = 2 is not supported anymore. If access to F is needed, then the material \n should overload giveFirstPKStressVector which has F as input.");
    }

    if ( this->nlGeometry != 0  &&  this->nlGeometry != 1 ) {
        OOFEM_ERROR("nlGeometry must be either 0 or 1 (%d not supported)", this->nlGeometry);
    } else {
        return 1;
    }
}

void
Rebar3d::NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    auto gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveCrossSection()->giveIPValue(answer, gp, type, tStep);
}


void
Rebar3d::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    FloatMatrix dN;
    this->interp.evaldNdx(dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 6);
    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 2) = dN.at(1, 2);
    answer.at(1, 3) = dN.at(1, 3);
    answer.at(1, 4) = dN.at(2, 1);
    answer.at(1, 5) = dN.at(2, 2);
    answer.at(1, 6) = dN.at(2, 3);
}


void
Rebar3d::computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    this->computeBmatrixAt(gp, answer);
}


void
Rebar3d::computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 1, this);
    }
}


double
Rebar3d::computeLength()
{
    return this->interp.giveLength(FEIElementGeometryWrapper(this) );
}


void
Rebar3d::computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // computes initial stress matrix of receiver (or geometric stiffness matrix)
    double area, NForce, contrib;
    double l0 = this->computeLength();
    FloatArray stress, strain;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    
    this->computeStrainVector(strain, gp, tStep);
    this->computeStressVector(stress, strain, gp, tStep);
    double Qc = this->giveCrossSectionReduction(gp, tStep, VM_Total);
    area = (1.0-Qc)*this->giveCrossSection()->give(CS_Area, gp);
    NForce = stress.at(1)*area;
    contrib = NForce / l0;
    
    answer.resize(6, 6);
    answer.zero();

    answer.at(2,2)=answer.at(3,3)=answer.at(5,5)=answer.at(6,6) = contrib;
    answer.at(2,5)=answer.at(5,2)=answer.at(3,6)=answer.at(6,3) = -contrib;
}

void
Rebar3d::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    answer.resize(6, 6);
    answer.zero();
    if ( !this->isActivated(tStep) ) {
        return;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double Qc = this->giveCrossSectionReduction(gp, tStep, VM_Total);
    double halfMass = density * (1.0-Qc)*this->giveCrossSection()->give(CS_Area, gp) * this->computeLength() * 0.5;
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(3, 3) = halfMass;
    answer.at(4, 4) = halfMass;
    answer.at(5, 5) = halfMass;
    answer.at(6, 6) = halfMass;
}


void
Rebar3d::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    this->interp.evalN(n, iLocCoord, FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 3);
}


double
Rebar3d::computeVolumeAround(GaussPoint *gp, TimeStep *tStep)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = this->interp.giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    double weight  = gp->giveWeight();
    double Qc = this->giveCrossSectionReduction(gp, tStep, VM_Total);
    return detJ * weight * (1.0-Qc)*this->giveCrossSection()->give(CS_Area, gp);
}

double
Rebar3d::computeMassCorroded(GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
// Returns the mass of the corroded part of the rebar element
{
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double Qc = this->giveCrossSectionReduction(gp, tStep, mode);
    double length = this->computeLength();
    double corrodedArea = Qc*this->giveCrossSection()->give(CS_Area, gp);
    return density * corrodedArea * length;
}

int
Rebar3d::giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx, ly(3), lz;

    lx.beDifferenceOf(this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    lx.normalize();

    ly(0) = lx(1);
    ly(1) = -lx(2);
    ly(2) = lx(0);

    // Construct orthogonal vector
    double npn = ly.dotProduct(lx);
    ly.add(-npn, lx);
    ly.normalize();
    lz.beVectorProductOf(ly, lx);

    answer.resize(3, 3);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


void
Rebar3d::computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveRealStress_1d(strain, gp, tStep);
}

int
Rebar3d::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_BeamForceMomentTensor ) {
        FloatArray stress, strain;
        double area;
        this->computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, gp, tStep);
        double Qc = this->giveCrossSectionReduction(gp, tStep, VM_Total);
        area = (1.0-Qc)*this->giveCrossSection()->give(CS_Area, gp);
        answer.resize(6);
        answer.at(1) = stress.at(1)*area;
        return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


void
Rebar3d::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_1d(rMode, gp, tStep);
}

void
Rebar3d::computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_dPdF_1d(rMode, gp, tStep);
}

void
Rebar3d::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_u, D_v, D_w };
}


void
Rebar3d::giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        OOFEM_ERROR("wrong edge number");
    }

    answer.resize(6);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
    answer.at(5) = 5;
    answer.at(6) = 6;
}

double
Rebar3d::computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        OOFEM_ERROR("wrong edge number");
    }

    double weight = gp->giveWeight();
    return this->interp.giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) * weight;
}


int
Rebar3d::computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatMatrix lcs;
    this->giveLocalCoordinateSystem(lcs);
    answer.beTranspositionOf(lcs);

    return 1;
}

double Rebar3d::giveCrossSectionReduction(GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
{
    double area_0 = this->giveCrossSection()->give(CS_Area, gp);
    double diameter_0 = 2.0 * sqrt(area_0 / M_PI);
    // get corrosion mass loss from the field
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr cf;
    if ( (cf = fm->giveField(FT_CorrosionMassLoss))) {
        FloatArray mloss, mloss1, mloss2;
        DofManager *dm1 = this->domain->giveDofManager(1);
        DofManager *dm2 = this->domain->giveDofManager(2);
        cf->evaluateAt(mloss1, dm1, mode, tStep);
        cf->evaluateAt(mloss2, dm2, mode, tStep);
        mloss= 0.5 * (mloss1 + mloss2); // average mass loss at the two nodes
        double density = this->giveCrossSection()->giveMaterial(gp)->give('d', gp);
        double x = mloss.at(1) / density; // corrosion mass loss to corrosion volume loss
        x = std::min(x, diameter_0/2.0); // clamp to max corrosion depth of radius
        double Q_c = 4.0 * (x / diameter_0- std::pow(x / diameter_0, 2.0));
        // Clamp Q_c to the range [0, 1]
        Q_c = std::min(std::max(Q_c, 0.0), 1.0);
        return Q_c;
        }
    return 0.0;
}    

#ifdef __OOFEG
void Rebar3d::drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p[ 2 ];  /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Rebar3d::drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p[ 2 ];  /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
