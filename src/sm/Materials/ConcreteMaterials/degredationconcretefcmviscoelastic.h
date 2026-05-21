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

#ifndef degredationconcretefcmviscoelastic_h
#define degredationconcretefcmviscoelastic_h

#include "degredationconcretefcm.h"
#include <memory>

#define _IFT_DegredationConcreteFCMViscoElastic_Name "degredationconcretefcmviscoelastic"
#define _IFT_ConcreteFCMViscoElastic_viscoMat "viscomat"
#define _IFT_ConcreteFCMViscoElastic_timedepfracturing "timedepfracturing"
#define _IFT_ConcreteFCMViscoElastic_fib_s "fib_s"
#define _IFT_ConcreteFCMViscoElastic_fcm28 "fcm28"
#define _IFT_ConcreteFCMViscoElastic_timeFactor "timefactor"
#define _IFT_ConcreteFCMViscoElastic_stiffnessFactor "stiffnessfactor"
#define _IFT_ConcreteFCMViscoElastic_gf28 "gf28"
#define _IFT_ConcreteFCMViscoElastic_ft28 "ft28"

namespace oofem {
class DegredationConcreteFCMViscoElasticStatus : public DegredationConcreteFCMStatus
{
protected:
    std :: unique_ptr< GaussPoint >slaveGpVisco;

    double var_ft = 0.;
    double var_gf = 0.;

public:
    DegredationConcreteFCMViscoElasticStatus(GaussPoint *g);

    double giveFractureEnergy() const { return var_gf; }
    void setFractureEnergy(double new_Gf) { var_gf = new_Gf; }

    double giveTensileStrength() const { return var_ft; }
    void setTensileStrength(double new_ft) { var_ft = new_ft; }

    GaussPoint *giveSlaveGaussPointVisco() { return this->slaveGpVisco.get(); }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "DegredationConcreteFCMViscoElasticStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


class DegredationConcreteFCMViscoElastic : public DegredationConcreteFCM
{
public:
    DegredationConcreteFCMViscoElastic(int n, Domain *d);

    void initializeFrom(const std::shared_ptr< InputRecord > &ir) override;
    const char *giveClassName() const override { return "DegredationConcreteFCMViscoElastic"; }
    const char *giveInputRecordName() const override { return _IFT_DegredationConcreteFCMViscoElastic_Name; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override { return std::make_unique<DegredationConcreteFCMViscoElasticStatus>(gp); }

    double give(int aProperty, GaussPoint *gp) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) const override;

    FloatArray computeStressIndependentStrainVector(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    MaterialStatus *giveStatus(GaussPoint *gp) const override;

protected:
    int viscoMat = 0;

    bool fib = false;
    double fib_s = 0.;
    double fib_fcm28 = 0.;
    double timeFactor = 0.;
    double stiffnessFactor = 0.;

    double giveTensileStrength(GaussPoint *gp, TimeStep *tStep) const override;
    double giveFractureEnergy(GaussPoint *gp, TimeStep *tStep) const override;

    double computeOverallElasticStiffness(GaussPoint *gp, TimeStep *tStep) const override;
    double computeOverallElasticShearModulus(GaussPoint *gp, TimeStep *tStep) const override;

    int checkConsistency(void) override;

    virtual double giveEquivalentTime(GaussPoint *gp, TimeStep *tStep) const;
};
} // end namespace oofem
#endif // degredationconcretefcmviscoelastic_h
