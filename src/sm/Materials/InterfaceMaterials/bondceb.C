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

#include "bondceb.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "fieldmanager.h"
#include "field.h"
#include "domain.h"
#include "engngm.h"


namespace oofem {
REGISTER_Material(BondCEBMaterial);

BondCEBMaterial :: BondCEBMaterial(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
{}


FloatArrayF<3>
BondCEBMaterial :: giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const
{
    BondCEBMaterialStatus *status = static_cast< BondCEBMaterialStatus * >( this->giveStatus(gp) );

    // normal traction evaluated elastically
    FloatArrayF<3> answer;
    answer.at(1) = kn * jump.at(1);

    // trial values of shear tractions evaluated elastically
    double s = 0., dKappa = 0.;
    for ( int i = 2; i <= 3; i++ ) {
        double depsi = jump.at(i) - status->giveJump().at(i);
        answer.at(i) = status->giveTraction().at(i) + ks * depsi;
        s += answer.at(i) * answer.at(i);
        dKappa += depsi*depsi;
     }

    // norm of trial shear traction
    s = sqrt(s);

    // cumulative slip increment
    dKappa = sqrt(dKappa);

    // cumulative slip at the end of the step
    double tempKappa = status->giveKappa() + dKappa;

    // maximum allowed norm of shear traction
    // get corrosion mass loss from the field
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr cf;
    double smax;
    if ( (cf = fm->giveField(FT_CorrosionFraction))) {
        double Qc = const_cast<BondCEBMaterial*>(this)->giveCrossSectionReduction(gp, tStep, VM_Total);
        double seq = 0.0136*Qc; // EC2020 13.6 mm*Qc with stirups
        double smax1 = evaluateBondStress(tempKappa); 
        double smax2 = evaluateBondStress(tempKappa+seq);
        smax = std::min(smax1, smax2);
        std::cout << "BondCEBMaterial" << " Qc: " << Qc << " smax1: " << smax1 << " smax2: " << smax2 << " smax: " << smax << std::endl;

    }else {
        smax = evaluateBondStress(tempKappa);
    }
    smax = evaluateBondStress(tempKappa);
    // reduce shear tractions, if needed
    if ( s > smax ) {
        for ( int i = 2; i <= 3; i++ ) {
            answer.at(i) *= smax / s;
        }
    }
 
    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
    status->setTempKappa(tempKappa);

    return answer;
}


FloatMatrixF<3,3>
BondCEBMaterial :: give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    ///@todo Only elastic tangent supported
    return diag<3>({kn, ks, ks});
}

double 
BondCEBMaterial :: giveCrossSectionReduction(GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
{
    double diameter_0 = diameter;
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
        double density = 7850; // density of steel in kg/m3
        double x = mloss.at(1) / density; // corrosion mass loss to corrosion volume loss
        double Q_c = 4.0 * (x / diameter_0- std::pow(x / diameter_0, 2.0));
        // Clamp Q_c to the range [0, 1]
        Q_c = std::min(std::max(Q_c, 0.0), 1.0);
        std::cout << "BondCEBMaterial" << " Q_c: " << Q_c << std::endl;
        std::cout << "BondCEBMaterial" << " mloss: " << mloss1.at(1) << " " << mloss2.at(1) << std::endl;
        return Q_c;
        }
    return 0.0;
}

double
BondCEBMaterial :: evaluateBondStress(const double kappa) const
{
    if ( kappa <= 0. )
        return 0.;
    if ( kappa <= s1 )
        return taumax * pow(kappa/s1, alpha);
    if ( kappa <= s2 )
        return taumax;
    if ( kappa <= s3 )
        return taumax - (taumax-tauf) * (kappa-s2) / (s3-s2);
    return tauf;
}


int
BondCEBMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    BondCEBMaterialStatus *status = static_cast< BondCEBMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveKappa();
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

void
BondCEBMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceMaterial :: initializeFrom(ir);

    // mandatory parameters
    IR_GIVE_FIELD(ir, kn, _IFT_BondCEBMaterial_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_BondCEBMaterial_ks);
    IR_GIVE_FIELD(ir, s1, _IFT_BondCEBMaterial_s1);
    IR_GIVE_FIELD(ir, s2, _IFT_BondCEBMaterial_s2);
    IR_GIVE_FIELD(ir, s3, _IFT_BondCEBMaterial_s3);
    IR_GIVE_FIELD(ir, taumax, _IFT_BondCEBMaterial_taumax);

    // optional parameters
    IR_GIVE_OPTIONAL_FIELD(ir, tauf, _IFT_BondCEBMaterial_tauf);
    IR_GIVE_OPTIONAL_FIELD(ir, alpha, _IFT_BondCEBMaterial_al);
    IR_GIVE_OPTIONAL_FIELD(ir, diameter, _IFT_BondCEBMaterial_diameter);

    // dependent parameter
    //@todo: why not simply if ks <taumax/s1 choose ks = taumax/s1
    //s0 = pow(pow(s1, -alpha)*taumax/ks, 1./(1.-alpha));
   if ( s0 > s1 ) {
     s0 = s1;
     ks = taumax/s1;
     OOFEM_WARNING("Parameter ks adjusted");
   }
}


void
BondCEBMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(this->kn, _IFT_BondCEBMaterial_kn);
    input.setField(this->ks, _IFT_BondCEBMaterial_ks);
    input.setField(this->s1, _IFT_BondCEBMaterial_s1);
    input.setField(this->s2, _IFT_BondCEBMaterial_s2);
    input.setField(this->s3, _IFT_BondCEBMaterial_s3);
    input.setField(this->taumax, _IFT_BondCEBMaterial_taumax);
}



BondCEBMaterialStatus :: BondCEBMaterialStatus(GaussPoint *g) : StructuralInterfaceMaterialStatus(g)
{}


void
BondCEBMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    StructuralInterfaceMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status { ");
    if ( this->kappa > 0.0 ) {
        fprintf(file, "kappa %g ", this->kappa);
    }

    fprintf(file, "}\n");
}


void
BondCEBMaterialStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();
    this->tempKappa = this->kappa;
}

void
BondCEBMaterialStatus :: updateYourself(TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
}


void
BondCEBMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: saveContext(stream, mode);

    if ( !stream.write(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

void
BondCEBMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    StructuralInterfaceMaterialStatus :: restoreContext(stream, mode);

    if ( !stream.read(kappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem
