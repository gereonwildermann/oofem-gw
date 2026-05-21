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

#ifndef rebar3d_h
#define rebar3d_h

#include "sm/Elements/structuralelement.h"
#include "sm/ErrorEstimators/directerrorindicatorrc.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"

#define _IFT_Rebar3d_Name "rebar3d"
#define _IFT_Rebar3d_nlgeoflag "nlgeo"

namespace oofem {
class FEI3dLineLin;
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;
class ParamKey;
/**
 * This class implements a two-node rebar element for three-dimensional
 * analysis.
 */
class Rebar3d : public StructuralElement,
    public ZZNodalRecoveryModelInterface,
    public NodalAveragingRecoveryModelInterface
{
protected:
    static FEI3dLineLin interp;
    int nlGeometry=0;
    static ParamKey IPK_Rebar3d_nlgeoflag;
public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Rebar3d(int n, Domain *d);
    virtual ~Rebar3d() { }

    /**
     * Returns the geometry mode describing the formulation used in the internal work
     * 0 - Engineering (small deformation) stress-strain mode
     * 1 - First Piola-Kirchhoff - Deformation gradient mode, P is defined as FS
     * 2 - Second Piola-Kirchhoff - Green-Lagrange strain mode with deformation gradient as input (deprecated and not supported)
     */
    int giveGeometryMode() { return nlGeometry; }
    /**
     * Computes the first Piola-Kirchhoff stress tensor on Voigt format. This method will
     * be called if nlGeo = 1 and mode = TL. This method computes the deformation gradient F and passes
     * it on to the crossection which then asks for the stress from the material.
     * @note P is related to S through F*S.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    void computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the Cauchy stress tensor on Voigt format. This method will
     * be called if nlGeo = 1 and mode = UL. This method computes the deformation gradient F and passes
     * it on to the crossection which then asks for the stress from the material.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    void computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the stiffness matrix of receiver.
     * The response is evaluated using @f$ \int B_{\mathrm{H}}^{\mathrm{T}} D B_{\mathrm{H}} \;\mathrm{d}v @f$, where
     * @f$ B_{\mathrm{H}} @f$ is the B-matrix which produces the displacement gradient vector @f$ H_{\mathrm{V}} @f$ when multiplied with
     * the solution vector a.
     * Reduced integration are taken into account.
     *
     * @param answer Computed stiffness matrix.
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;



    /**
     * Computes the initial stiffness matrix of receiver. This method is used only if mode = UL
     * The response is evaluated using @f$ \int B ({\mathrm{\sigma}}\otimes \delta )B_{\mathrm{H}} \;\mathrm{d}v @f$, where
     * @f$ B @f$ is the classical B-matrix, but computed wrt updated node position
     *
     * @param answer Computed initial stiffness matrix.
     * @param tStep Time step.
     */
    void computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep) override;

    /**
     * Computes the stiffness matrix of receiver.
     * The response is evaluated using @f$ \int B_{\mathrm{H}}^{\mathrm{T}} D B_{\mathrm{H}} \;\mathrm{d}v @f$, where
     * @f$ B_{\mathrm{H}} @f$ is the B-matrix which produces the displacement gradient vector @f$ H_{\mathrm{V}} @f$ when multiplied with
     * the solution vector a.
     * @note Reduced intergration is not taken into account.
     * The integration procedure uses an integrationRulesArray for numerical integration. Each integration rule is
     * considered to represent a separate sub-cell/element. Typically this would be used when integration of the element
     * domain needs special treatment, e.g. when using the XFEM.
     *
     * @param answer Computed stiffness matrix.
     * @param rMode Response mode.
     * @param tStep Time step.
     */
    void computeStiffnessMatrix_withIRulesAsSubcells(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    /**
     * Evaluates nodal representation of real internal forces.
     * Necessary transformations are taken into account. @todo what is meant?
     *
     * @param answer Equivalent nodal forces vector.
     * @param tStep Time step
     * @param useUpdatedGpRecord If equal to zero, the stresses in integration points are computed (slow but safe).
     */

    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    /**
     * Computes large strain constitutive matrix of receiver. Default implementation uses element cross section
     * giveCharMaterialStiffnessMatrix service.
     * @param answer Constitutive matrix.
     * @param rMode Material response mode of answer.
     * @param gp Integration point for which constitutive matrix is computed.
     * @param tStep Time step.
     */
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer,
                                                   MatResponseMode rMode, GaussPoint *gp,
                                                   TimeStep *tStep);


    /**
     * Evaluates nodal representation of real internal forces.
     *
     * Numerical integration procedure uses integrationRulesArray
     * for numerical integration. The integration procedure uses an integrationRulesArray for numerical integration.
     * Each integration rule is considered to represent a separate sub-cell/element. Typically this would be used when
     * integration of the element domain needs special treatment, e.g. when using the XFEM.
     *
     * @param answer Equivalent nodal forces vector.
     * @param tStep Time step.
     * @param useUpdatedGpRecord If equal to zero, the stresses in the integration points are computed (slow but safe).
     */
    void giveInternalForcesVector_withIRulesAsSubcells(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;

    /**
     * Computes the deformation gradient in Voigt form at integration point ip and at time
     * step tStep. Computes the displacement gradient and adds an identitiy tensor.
     *
     * @param answer Deformation gradient vector
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    virtual void computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the current volume of element
     */
    double computeCurrentVolume(TimeStep *tStep);
    double giveCrossSectionReduction(GaussPoint *gp, TimeStep *tStep, ValueModeType mode = VM_Total);

    FEInterpolation *giveInterpolation() const override;

    double computeLength() override;

    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override
    { this->computeLumpedMassMatrix(answer, tStep); }
    int giveLocalCoordinateSystem(FloatMatrix &answer) override;

    int computeNumberOfDofs() override { return 6; }
    void giveDofManDofIDMask(int inode, IntArray &) const override;


    // characteristic length (for crack band approach)
    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override
    { return this->computeLength(); }

    double computeVolumeAround(GaussPoint *gp, TimeStep *tStep);
    double computeMassCorroded(GaussPoint *gp, TimeStep *tStep, ValueModeType mode = VM_Total);

    int testElementExtension(ElementExtension ext) override { return ( ext == Element_EdgeLoadSupport ); }

    Interface *giveInterface(InterfaceType it) override;

    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep) override;

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) override;
#endif

    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Rebar3d_Name; }
    const char *giveClassName() const override { return "Rebar3d"; }
    void initializeFrom(const std::shared_ptr<InputRecord> &ir, int priority) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    MaterialMode giveMaterialMode() override { return _1dMat; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_line_1;}

    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

protected:
    int checkConsistency() override;
    /**
     * Computes a matrix which, multiplied by the column matrix of nodal displacements,
     * gives the displacement gradient stored by columns.
     * The components of this matrix are derivatives of the shape functions,
     * but they are arranged in a somewhat different way from the usual B matrix.
     * @param gp Integration point.
     * @param answer BF matrix at this point.
     */
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    // edge load support
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    double computeEdgeVolumeAround(GaussPoint *gp, int) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *gp) override;
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS) override;
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // rebar3d_h
