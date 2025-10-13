/*
 *  OOFEM : Object Oriented Finite Element Code
 *  RebarCrossSection - circular cross-section with optional corrosion loss
 *
 *  Copyright (C) 2025  Gereon Wildermann
 *
 *  Licensed under the GNU Lesser General Public License v2.1 or later.
 */

#include "sm/CrossSections/rebarcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Elements/nlstructuralelement.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "core/field.h"
#include "fieldmanager.h"
#include "material.h"
#include "engngm.h"
#include "domain.h"

#include <cmath>

namespace oofem {
REGISTER_CrossSection(RebarCrossSection);

void RebarCrossSection::initializeFrom(InputRecord &ir)
{
    // read diameter if present
    if ( ir.hasField(_IFT_RebarCrossSection_diameter) ) {
        IR_GIVE_FIELD(ir, diameter, _IFT_RebarCrossSection_diameter);
    }

    // call parent for any other initialization
    SimpleCrossSection::initializeFrom(ir);
}

void RebarCrossSection::giveInputRecord(DynamicInputRecord &input)
{
    // call parent first
    SimpleCrossSection::giveInputRecord(input);

    // add our diameter field
    input.setField(this->diameter, _IFT_RebarCrossSection_diameter);
}

double RebarCrossSection::give(CrossSectionProperty aProperty, GaussPoint *gp) const
{
    if (aProperty == CS_Area) {
        double A0 = M_PI * diameter * diameter / 4.0;

        if (!gp) // no Gauss point, return full area
            return A0;

        auto mat = dynamic_cast< StructuralMaterial * >( this->giveMaterial(gp) );
        if (!mat)
            return A0;

        // get corrosion mass loss from the field
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
        FieldPtr cf;
        if ( ( cf = fm->giveField(FT_CorrosionMassLoss)) ) {
            FloatArray gcoords, mloss;
            int err;
            Element *elem = gp->giveElement()
            elem ->computeGlobalCoordinates(gcoords, gp->giveNaturalCoordinates());
            if (( err = cf->evaluateAt(mloss, gcoords, mode, tStep))) {
                OOFEM_ERROR("cf->evaluateAt failed, element %d, error code %d", gp->giveElement()->giveNumber(), err);
            }else {
                if (mLoss <= 0.0)
                return A0;

            double rho = mat->give('d', gp);    // material density
            double x = mLoss / rho;                  // corrosion depth
            double Qc = 4.0 * (x / diameter - std::pow(x / diameter, 2.0));
            Qc = std::min(std::max(Qc, 0.0), 1.0);
            return (1.0 - Qc) * A0;
        } else {
            return A0;
        }

    // default: call parent
    return SimpleCrossSection::give(aProperty, gp);
}

} // namespace oofem
