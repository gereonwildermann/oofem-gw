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
 *               Copyright (C) 1993 - 2016   Borek Patzak
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

#ifndef degredationconcretefcm_h
#define degredationconcretefcm_h

#include "sm/Materials/fcm.h"
#include "randommaterialext.h"

///@name Input fields for ConcreteFCM
//@{
#define _IFT_DegredationConcreteFCM_Name "degredationconcretefcm"
#define _IFT_ConcreteFCM_softType "softtype"
#define _IFT_ConcreteFCM_shearType "sheartype"
#define _IFT_ConcreteFCM_shearStrengthType "shearstrengthtype"
#define _IFT_ConcreteFCM_gf "gf"
#define _IFT_ConcreteFCM_ft "ft"
#define _IFT_ConcreteFCM_beta "beta"
#define _IFT_ConcreteFCM_sf "sf"
#define _IFT_ConcreteFCM_sf_numer "sf_numer"
#define _IFT_ConcreteFCM_fc "fc"
#define _IFT_ConcreteFCM_ag "ag"
#define _IFT_ConcreteFCM_lengthScale "lengthscale"
#define _IFT_ConcreteFCM_soft_w "soft_w"
#define _IFT_ConcreteFCM_soft_function_w "soft(w)"
#define _IFT_ConcreteFCM_soft_eps "soft_eps"
#define _IFT_ConcreteFCM_soft_function_eps "soft(eps)"
#define _IFT_ConcreteFCM_beta_w "beta_w"
#define _IFT_ConcreteFCM_beta_function "beta(w)"
#define _IFT_ConcreteFCM_H "h"
#define _IFT_ConcreteFCM_eps_f "eps_f"
#define _IFT_ConcreteFCM_beta_Gf "beta_Gf"
#define _IFT_ConcreteFCM_beta_f "beta_f"
#define _IFT_ConcreteFCM_beta_E "beta_E"
//@}

namespace oofem {
class DegredationConcreteFCMStatus : public FCMMaterialStatus, public RandomMaterialStatusExtensionInterface
{
protected:
    double elasticStiffnessScale = 1.;
    double tempElasticStiffnessScale = 1.;
    int lastStressRescalingStep = -1;

public:
    DegredationConcreteFCMStatus(GaussPoint *g);

    double giveElasticStiffnessScale() const { return elasticStiffnessScale; }
    double giveTempElasticStiffnessScale() const { return tempElasticStiffnessScale; }
    void setTempElasticStiffnessScale(double value) { tempElasticStiffnessScale = value; }

    int giveLastStressRescalingStep() const { return lastStressRescalingStep; }
    void setLastStressRescalingStep(int step) { lastStressRescalingStep = step; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "DegredationConcreteFCMStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    Interface *giveInterface(InterfaceType it) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


class DegredationConcreteFCM : public FCMMaterial, public RandomMaterialExtensionInterface
{
public:
    DegredationConcreteFCM(int n, Domain *d);

    void initializeFrom(const std::shared_ptr< InputRecord > &ir) override;
    const char *giveClassName() const override { return "DegredationConcreteFCM"; }
    const char *giveInputRecordName() const override { return _IFT_DegredationConcreteFCM_Name; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<DegredationConcreteFCMStatus>(gp); }

    double give(int aProperty, GaussPoint *gp) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) const override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
    double Gf = 0.;
    double Ft = 0.;
    double beta = 0.;
    double sf = 0.;
    double sf_numer = 0.;

    double beta_Gf = 1.;
    double beta_f = 1.;
    double beta_E = 1.;

    double fc = 0.;
    double ag = 0.;
    double lengthScale = 0.;

    FloatArray soft_w, soft_function_w;
    FloatArray soft_eps, soft_function_eps;
    FloatArray beta_w, beta_function;

    double H = 0.;
    double eps_f = 0.;

    double giveTensileStrength(GaussPoint *gp, TimeStep *tStep) const override;
    virtual double giveFractureEnergy(GaussPoint *gp, TimeStep *tStep) const;
    virtual double computeOverallElasticStiffness(GaussPoint *gp, TimeStep *tStep) const override;
    // virtual double computeOverallElasticShearModulus(GaussPoint *gp, TimeStep *tStep) const override;

    double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, int i) const override;
    double giveCrackingModulusInTension(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, int i) const override;
    double computeEffectiveShearModulus(GaussPoint *gp, TimeStep *tStep, int i) const override;
    double computeD2ModulusForCrack(GaussPoint *gp, TimeStep *tStep, int icrack) const override;
    double computeNumerD2ModulusForCrack(GaussPoint *gp, TimeStep *tStep, int icrack) const override;
    double giveNormalCrackingStress(GaussPoint *gp, TimeStep *tStep, double eps_cr, int i) const override;
    double maxShearStress(GaussPoint *gp, TimeStep *tStep, int i) const override;

    virtual double computeResidualTensileStrength(GaussPoint *gp, TimeStep *tStep) const;

    void checkSnapBack(GaussPoint *gp, TimeStep *tStep, int crack) const override;

    enum SofteningType { ST_NONE, ST_Exponential, ST_Linear, ST_Hordijk, ST_UserDefinedCrack, ST_LinearHardeningStrain, ST_UserDefinedStrain, ST_Unknown };
    SofteningType softType = ST_NONE;

    enum ShearRetentionType { SHR_NONE, SHR_Const_ShearRetFactor, SHR_Const_ShearFactorCoeff, SHR_UserDefined_ShearRetFactor, SHR_Unknown };
    ShearRetentionType shearType = SHR_Unknown;

    enum ShearStrengthType { SHS_NONE, SHS_Const_Ft, SHS_Collins_Interlock, SHS_Residual_Ft, SHS_Unknown };
    ShearStrengthType shearStrengthType = SHS_Unknown;
};
} // end namespace oofem
#endif // degredationconcretefcm_h
