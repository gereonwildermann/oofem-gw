/*
 *  OOFEM : Object Oriented Finite Element Code
 *  RebarCrossSection - circular cross-section with optional corrosion loss
 *
 *  Copyright (C) 2025  Gereon Wildermann
 *
 *  Licensed under the GNU Lesser General Public License v2.1 or later.
 */

#ifndef rebarcrosssection_h
#define rebarcrosssection_h

#include "sm/CrossSections/simplecrosssection.h"
#include "floatarray.h"

/// @name Input fields for RebarCrossSection
//@{
#define _IFT_RebarCrossSection_Name "rebarcs"
#define _IFT_RebarCrossSection_diameter "diameter"
//@}

namespace oofem {
/**
 * RebarCrossSection:
 * - If corrosion mass loss is not defined, computes A = π d² / 4.
 * - If corrosion mass loss (field FT_CorrosionMassLoss) is defined:
 *     * get density ρ from the associated material,
 *     * compute corrosion depth x = m_loss / ρ,
 *     * compute area reduction factor Q_c = 4 (x/d0 - (x/d0)^2),
 *     * compute reduced area A = (1 - Q_c) A0.
 */
class OOFEM_EXPORT RebarCrossSection : public SimpleCrossSection
{
protected:
    double diameter;  ///< Initial bar diameter

public:
    /// Constructor
    RebarCrossSection(int n, Domain *d) : SimpleCrossSection(n, d), diameter(0.0) {}

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    const char *giveClassName() const override { return "RebarCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_RebarCrossSection_Name; }

    double give(CrossSectionProperty a, GaussPoint *gp) const override;
};
} // end namespace oofem

#endif // rebarcrosssection_h
