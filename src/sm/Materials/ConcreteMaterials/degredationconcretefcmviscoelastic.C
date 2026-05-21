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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#include "degredationconcretefcmviscoelastic.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include "contextioerr.h"
#include "../RheoChainMaterials/rheoChM.h"
#include "fieldmanager.h"
#include "field.h"
#include "domain.h"
#include "engngm.h"
#include "sm/Elements/structuralelement.h"

namespace {
double evaluateDegreeOfDegradation(oofem::Domain *domain, oofem::GaussPoint *gp, oofem::TimeStep *tStep)
{
    oofem::FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    oofem::FieldPtr ddf;
    if ( !( ddf = fm->giveField(oofem::FT_DegreeOfDegradation) ) ) {
        return 0.;
    }

    oofem::Coordinates gcoords, lcoords;
    oofem::FloatArray dod;
    lcoords = gp->giveNaturalCoordinates();
    int err = 0;
    static_cast< oofem::StructuralElement * >( gp->giveElement() )->computeGlobalCoordinates(gcoords, lcoords);
    if ( ( err = ddf->evaluateAt(dod, gcoords, oofem::VM_Total, tStep) ) ) {
        throw oofem::RuntimeException(__func__, __FILE__, __LINE__, "ddf->evaluateAt failed, element %d, error code %d", gp->giveElement()->giveNumber(), err);
    }

    if ( dod.giveSize() < 1 ) {
        return 0.;
    }

    return oofem::min(1., oofem::max(0., dod.at(1)));
}
}

namespace oofem {
REGISTER_Material(DegredationConcreteFCMViscoElastic);

DegredationConcreteFCMViscoElastic::DegredationConcreteFCMViscoElastic(int n, Domain *d) : DegredationConcreteFCM(n, d)
{}


void
DegredationConcreteFCMViscoElastic::initializeFrom(const std::shared_ptr< InputRecord > &ir)
{
    DegredationConcreteFCM::initializeFrom(ir);

    IR_GIVE_FIELD(ir, viscoMat, _IFT_ConcreteFCMViscoElastic_viscoMat);

    this->fib = false;

    if ( ir->hasField(_IFT_ConcreteFCMViscoElastic_timedepfracturing) ) {
        this->fib = true;
        //
        IR_GIVE_FIELD(ir, fib_s, _IFT_ConcreteFCMViscoElastic_fib_s);
        // the same compressive strength as for the prediction using the B3 formulas
        IR_GIVE_FIELD(ir, fib_fcm28, _IFT_ConcreteFCMViscoElastic_fcm28);

        IR_GIVE_FIELD(ir, timeFactor, _IFT_ConcreteFCMViscoElastic_timeFactor);
        IR_GIVE_FIELD(ir, stiffnessFactor, _IFT_ConcreteFCMViscoElastic_stiffnessFactor);

        DegredationConcreteFCM::Gf = -1.;
        DegredationConcreteFCM::Ft = -1.;
        IR_GIVE_OPTIONAL_FIELD(ir, DegredationConcreteFCM::Ft, _IFT_ConcreteFCMViscoElastic_ft28);
        IR_GIVE_OPTIONAL_FIELD(ir, DegredationConcreteFCM::Gf, _IFT_ConcreteFCMViscoElastic_gf28);
    }


    if ( propertyDictionary.includes(tAlpha) ) {
        if ( propertyDictionary.at(tAlpha) != 0. ) {
            OOFEM_ERROR("tAlpha must be set to zero in ConcreteFCMViscoElastic material");
        }
    }
}


double
DegredationConcreteFCMViscoElastic::give(int aProperty, GaussPoint *gp) const
{
    return DegredationConcreteFCM::give(aProperty, gp);
}


void
DegredationConcreteFCMViscoElastic::giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                              const FloatArray &totalStrain,
                                              TimeStep *tStep) const
{
    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    FloatArray partialStrain;
    FloatArray viscoStress;

    FCMMaterial::giveRealStressVector(answer, gp, totalStrain, tStep);

    partialStrain = totalStrain;

    if ( status->giveNumberOfTempCracks() > 0 ) {
        FloatArray crackStrain;
        FloatMatrix epsL2G = status->giveL2GStrainVectorTransformationMtrx();
        // from local to global
        //    crackStrain.rotatedWith(epsL2G, 'n');
        crackStrain.beProductOf( epsL2G, status->giveTempCrackStrainVector() );
        partialStrain.subtract(crackStrain);
    }

    // viscoStress should be equal to answer - this method is called only for updating
    rheoMat->giveRealStressVector(viscoStress, status->giveSlaveGaussPointVisco(), partialStrain, tStep);

    return;
}

FloatArray
DegredationConcreteFCMViscoElastic::computeStressIndependentStrainVector(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
{
    // temperature strain is treated ONLY by the rheoMat

    if ( !this->isActivated(tStep) ) {
        return FloatArray();
    }

    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    return rheoMat->computeStressIndependentStrainVector(status->giveSlaveGaussPointVisco(), tStep, mode);
}


int
DegredationConcreteFCMViscoElastic::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    if ( type == IST_TensileStrength ) {
        answer.resize(1);
        answer.at(1) = status->giveTensileStrength();
        return 1;
    } else if ( type == IST_ResidualTensileStrength ) {
        double sigma;
        int nCracks;
        double emax;

        nCracks = status->giveNumberOfTempCracks();
        sigma = status->giveTensileStrength();

        for ( int i = 1; i <= nCracks; i++ ) {
            emax = status->giveMaxCrackStrain(i);

            if ( emax > 0. ) {
                // normalisation wrt crack number is done inside the function
                // emax /= this->giveNumberOfCracksInDirection(gp, i);

                sigma = this->giveNormalCrackingStress(gp, tStep, emax, i);
                sigma = min(sigma, this->giveTensileStrength(gp, tStep) );
            }
        }

        answer.resize(1);
        answer.zero();
        answer.at(1) =  sigma;

        return 1;
    }  else if ( type == IST_CrackIndex ) {
        answer.resize(1);
        answer.zero();

        // cracking is initiated
        if ( status->giveNumberOfCracks() ) {
            answer.at(1) = 1.;
        } else {
            FloatArray sigma;
            FloatArray princStress;
            StructuralMaterial::giveFullSymVectorForm( sigma, status->giveStressVector(), gp->giveMaterialMode() );
            this->computePrincipalValues(princStress, sigma, principal_stress);
            answer.at(1) = max( 0., princStress.at(1) / status->giveTensileStrength() );
        }
        return 1;
    }

    //   if ( ( type == IST_DryingShrinkageTensor ) ||
    //        ( type == IST_AutogenousShrinkageTensor ) ||
    //        ( type == IST_TotalShrinkageTensor ) ||
    //        ( type == IST_CreepStrainTensor ) ||
    //        ( type == IST_DryingShrinkageTensor ) ||
    //        ( type == IST_ThermalStrainTensor ) ) {

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
    if ( rheoMat->giveIPValue(answer, status->giveSlaveGaussPointVisco(), type, tStep) ) {
        return 1;
    }

    return DegredationConcreteFCM::giveIPValue(answer, gp, type, tStep);
}



MaterialStatus* 
DegredationConcreteFCMViscoElastic::giveStatus(GaussPoint *gp) const
{
    if (gp->hasMaterialStatus()) {
        return static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    } else {
        MaterialStatus *status = static_cast<MaterialStatus*> (gp->setMaterialStatus (this->CreateStatus(gp)));
        this->_generateStatusVariables(gp);
        return status;
    }
}


double
DegredationConcreteFCMViscoElastic::giveTensileStrength(GaussPoint *gp, TimeStep *tStep) const
{
    // return value
    double ftm = 0.;
    const double DoD = evaluateDegreeOfDegradation(domain, gp, tStep);

    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    // check if the fracture properties are constant or time-dependent
    if ( fib ) {
        // compressive strength
        double fcm, fcm28mod;
        double equivalentTime;

        // virgin material-> compute!
        if ( status->giveNumberOfTempCracks() == 0 ) {
            equivalentTime = this->giveEquivalentTime(gp, tStep);

            if ( DegredationConcreteFCM::Ft > 0. ) { // user-specified 28-day strength
                fcm28mod = pow(DegredationConcreteFCM::Ft * this->stiffnessFactor / 0.3e6, 3. / 2.) + 8.;
                fcm = exp(fib_s * ( 1. - sqrt(28. * this->timeFactor / equivalentTime) ) ) * fcm28mod;
            } else { // returns fcm in MPa - formula 5.5-51
                fcm = exp(this->fib_s * ( 1. - sqrt(28. * this->timeFactor / equivalentTime) ) ) * this->fib_fcm28;
            }

            // ftm adjusted according to the stiffnessFactor (MPa by default)
            if ( fcm >= 58. ) {
                ftm = 2.12 * log(1. + 0.1 * fcm) * 1.e6 / this->stiffnessFactor;
            } else if ( fcm <= 20. ) {
                ftm = 0.07862 * fcm * 1.e6 / this->stiffnessFactor; // 12^(2/3) * 0.3 / 20 = 0.07862
            } else {
                ftm = 0.3 * pow(fcm - 8., 2. / 3.) * 1.e6 / this->stiffnessFactor; //5.1-3a
            }

            // ftm adjusted according to the stiffnessFactor (MPa by default)
            /*      if ( fcm >= 20. ) {
             * ftm = 0.3 * pow(fcm - 8., 2. / 3.) * 1.e6 / this->stiffnessFactor; //5.1-3a
             * } else if ( fcm < 8. ) {
             * // upper formula does not hold for concretes with fcm < 8 MPa
             * ftm = 0.3 * pow(fcm, 2. / 3.) * 1.e6 / this->stiffnessFactor;
             * } else {
             * // smooth transition
             * ftm = 0.3 * pow(fcm - ( 8. * ( fcm - 8. ) / ( 20. - 8. ) ), 2. / 3.) * 1.e6 / this->stiffnessFactor;
             * }*/

            status->setTensileStrength(ftm);
        } else {
            ftm = status->giveTensileStrength();
        }
    } else {
        ftm = DegredationConcreteFCM::giveTensileStrength(gp, tStep);
    }

    ftm *= ( 1. - ( 1. - beta_f ) * DoD );
    status->setTensileStrength(ftm);

    return ftm;
}

double
DegredationConcreteFCMViscoElastic::giveFractureEnergy(GaussPoint *gp, TimeStep *tStep) const
{
    // return value
    double Gf = 0.;
    const double DoD = evaluateDegreeOfDegradation(domain, gp, tStep);
    double Gf28;
    double ftm, ftm28;

    // check if the fracture properties are constant or time-dependent
    if ( fib ) {
        DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

        // virgin material-> compute!
        if ( status->giveNumberOfTempCracks() == 0 ) {
            // 1)
            if ( DegredationConcreteFCM::Gf > 0. ) {
                Gf28 = DegredationConcreteFCM::Gf;
            } else {
                Gf28 = 73. * pow(fib_fcm28, 0.18) / this->stiffnessFactor;
            }

            // 2)
            ftm = this->giveTensileStrength(gp, tStep);

            if ( DegredationConcreteFCM::Ft > 0. ) { // user-specified 28-day strength
                ftm28 = DegredationConcreteFCM::Ft;
            } else {
                if ( fib_fcm28 >= 58. ) {
                    ftm28 = 2.12 * log(1. + 0.1 * fib_fcm28) * 1.e6 / this->stiffnessFactor;
                } else if ( fib_fcm28 <= 20. ) {
                    ftm28 = 0.07862 * fib_fcm28 * 1.e6 / this->stiffnessFactor; // 12^(2/3) * 0.3 / 20 = 0.07862
                } else {
                    ftm28 = 0.3 * pow(fib_fcm28 - 8., 2. / 3.) * 1.e6 / this->stiffnessFactor; //5.1-3a
                }
            }

            // 3)
            Gf = Gf28 * ftm / ftm28;

            status->setFractureEnergy(Gf);
        } else {
            Gf = status->giveFractureEnergy();
        }
    } else {
        Gf = DegredationConcreteFCM::giveFractureEnergy(gp, tStep);
    }

    Gf *= ( 1. - ( 1. - beta_Gf ) * DoD );

    return Gf;
}


double
DegredationConcreteFCMViscoElastic::computeOverallElasticStiffness(GaussPoint *gp, TimeStep *tStep) const {
    double stiffness;
    const double DoD = evaluateDegreeOfDegradation(domain, gp, tStep);

    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    stiffness = rheoMat->giveEModulus(status->giveSlaveGaussPointVisco(), tStep) * ( 1. - ( 1. - beta_E ) * DoD );

    return stiffness;
}

double
DegredationConcreteFCMViscoElastic::computeOverallElasticShearModulus(GaussPoint *gp, TimeStep *tStep) const {
    double Evisco;
    double Gvisco;

    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    Evisco = rheoMat->giveEModulus(status->giveSlaveGaussPointVisco(), tStep);

    Gvisco = FCMMaterial::linearElasticMaterial.giveShearModulus() * Evisco / this->give('E', gp);

    return Gvisco;
}



int DegredationConcreteFCMViscoElastic::checkConsistency()
{
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );

    MaterialMode helpMM;
    helpMM = MaterialMode(_Unknown);
    GaussPoint helpGP(NULL, 0, 0., helpMM);


    if ( rheoMat->givePoissonsRatio() != this->linearElasticMaterial.givePoissonsRatio() ) {
        OOFEM_ERROR("The Poisson ratio of the fracturing and viscoelastic material are not equal.");
    }

    if ( fabs(this->linearElasticMaterial.give(tAlpha, & helpGP) ) != 0. ) {
        OOFEM_ERROR("tAlpha must be set to zero in ConcreteFCMViscoElastic material");
    }

    return FEMComponent::checkConsistency();
}

double
DegredationConcreteFCMViscoElastic::giveEquivalentTime(GaussPoint *gp, TimeStep *tStep) const
{
    RheoChainMaterial *rheoMat = static_cast< RheoChainMaterial * >( domain->giveMaterial(this->viscoMat) );
    DegredationConcreteFCMViscoElasticStatus *status = static_cast< DegredationConcreteFCMViscoElasticStatus * >( gp->giveMaterialStatus() );

    return rheoMat->giveEquivalentTime(status->giveSlaveGaussPointVisco(), tStep);
}




///////////////////////////////////////////////////////////////////
//                      CONCRETE FCM STATUS                     ///
///////////////////////////////////////////////////////////////////

DegredationConcreteFCMViscoElasticStatus::DegredationConcreteFCMViscoElasticStatus(GaussPoint *gp) :
    DegredationConcreteFCMStatus(gp),
    slaveGpVisco( std::make_unique< GaussPoint >(gp->giveIntegrationRule(), 0, gp->giveNaturalCoordinates(), 0., gp->giveMaterialMode() ) )
{}


void
DegredationConcreteFCMViscoElasticStatus::printOutputAt(FILE *file, TimeStep *tStep) const
{
    DegredationConcreteFCMStatus::printOutputAt(file, tStep);

    fprintf(file, "remark {Output for slave viscoelastic material}\n");
    this->slaveGpVisco->giveMaterialStatus()->printOutputAt(file, tStep);
    fprintf(file, "\n");
}


void
DegredationConcreteFCMViscoElasticStatus::initTempStatus()
//
// initializes temp variables according to variables form previous equlibrium state.
//
{
    DegredationConcreteFCMStatus::initTempStatus();

    RheoChainMaterialStatus *rheoStatus = static_cast< RheoChainMaterialStatus * >( this->giveSlaveGaussPointVisco()->giveMaterialStatus() );

    rheoStatus->initTempStatus();
}



void
DegredationConcreteFCMViscoElasticStatus::updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables correspond to newly reched equilibrium.
//
{
    DegredationConcreteFCMStatus::updateYourself(tStep);

    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->updateYourself(tStep);
}


void
DegredationConcreteFCMViscoElasticStatus::saveContext(DataStream &stream, ContextMode mode)
//
// saves full information stored in this Status
// no temp variables stored
//
{
    // save cracking status
    DegredationConcreteFCMStatus::saveContext(stream, mode);
    // save viscoelastic status
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->saveContext(stream, mode);
}

void
DegredationConcreteFCMViscoElasticStatus::restoreContext(DataStream &stream, ContextMode mode)
//
// restores full information stored in stream to this Status
//
{
    // read parent cracking class status
    DegredationConcreteFCMStatus::restoreContext(stream, mode);

    // restore viscoelastic amterial
    this->giveSlaveGaussPointVisco()->giveMaterialStatus()->restoreContext(stream, mode);
}
} // end namespace oofem
