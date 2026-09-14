//==============================================================================
/*
    Tests for the iterative refinement of the OBB axes.
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBBFit.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
#include <random>
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

static double product(const cVector3d& a_size)
{
    return (a_size(0) * a_size(1) * a_size(2));
}

static double axisAngleDeg(const cVector3d& a_dir, const cVector3d& a_axis)
{
    double c = cClamp(fabs(cDot(cNormalize(a_dir), cNormalize(a_axis))), 0.0, 1.0);
    return (cRadToDeg(acos(c)));
}

static cMatrix3d rotationDeg(const cVector3d& a_axis, double a_angleDeg)
{
    cMatrix3d rot;
    rot.setAxisAngleRotationDeg(a_axis, a_angleDeg);
    return (rot);
}

static std::vector<cVector3d> boxCorners(const cVector3d& a_half, const cMatrix3d& a_rot, const cVector3d& a_pos)
{
    std::vector<cVector3d> pts;
    for (int i=0; i<8; i++)
        pts.push_back(a_rot * cVector3d((i & 1) ? a_half(0) : -a_half(0),
                                        (i & 2) ? a_half(1) : -a_half(1),
                                        (i & 4) ? a_half(2) : -a_half(2)) + a_pos);
    return (pts);
}

static bool containsAll(const cCollisionOBBBox& a_box, const std::vector<cVector3d>& a_points)
{
    for (size_t i=0; i<a_points.size(); i++)
        if (!a_box.contains(a_points[i])) return (false);
    return (true);
}

static cOBBRefineSettings noRefinement()
{
    cOBBRefineSettings settings;
    settings.m_enabled = false;
    return (settings);
}


int main()
{
    std::mt19937 rng(7);

    //--------------------------------------------------------------------------
    // box size for given axes
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> pts = boxCorners(cVector3d(2, 1, 0.5), cIdentity3d(), cVector3d(3, 3, 3));
        CHECK(approx(cComputeBoxSize(pts, cIdentity3d()), cVector3d(4, 2, 1), 1e-12));
        cVector3d size45 = cComputeBoxSize(pts, rotationDeg(cVector3d(0, 0, 1), 45));
        // CHAI3D's degree-to-radian constant has ~11 significant digits
        CHECK(approx(size45, cVector3d(3 * sqrt(2.0), 3 * sqrt(2.0), 1), 1e-9));
    }

    //--------------------------------------------------------------------------
    // recovery of the true axes from a perturbed frame
    //--------------------------------------------------------------------------
    {
        const cMatrix3d R = rotationDeg(cVector3d(1, 2, 3), 37);
        std::vector<cVector3d> pts = boxCorners(cVector3d(2, 1, 0.5), R, cVector3d(10, -5, 3));

        // frame off by 3, -2 and 4 degrees about its own axes
        cMatrix3d perturbed = R * rotationDeg(cVector3d(1, 0, 0), 3)
                                * rotationDeg(cVector3d(0, 1, 0), -2)
                                * rotationDeg(cVector3d(0, 0, 1), 4);

        cOBBRefineStats stats;
        cMatrix3d refined = cRefineOBBAxes(pts, perturbed, cOBBRefineSettings(), &stats);
        double err = cMax(axisAngleDeg(refined.getCol0(), R.getCol0()),
                     cMax(axisAngleDeg(refined.getCol1(), R.getCol1()),
                          axisAngleDeg(refined.getCol2(), R.getCol2())));

        CHECK(approx(stats.m_initialVolume, product(cComputeBoxSize(pts, perturbed)), 1e-9));
        CHECK(approx(stats.m_finalVolume, product(cComputeBoxSize(pts, refined)), 1e-9));
        // with a finest step of 1 degree the residual error is below 1 degree
        CHECK(stats.m_finalVolume < stats.m_initialVolume);
        CHECK(stats.m_finalVolume < 8.0 * 1.05);
        CHECK(err < 1.0);
        printf("perturbed frame, steps {5,1} deg:     volume %.4f -> %.4f (exact 8), axis error %.3f deg, "
               "%d iterations, %d rotations, %d evaluations\n",
               stats.m_initialVolume, stats.m_finalVolume, err,
               stats.m_iterations, stats.m_acceptedRotations, stats.m_evaluations);

        // a finer last step gives a finer result
        cOBBRefineSettings fine;
        fine.m_anglesDeg.push_back(0.2);
        cOBBRefineStats fineStats;
        cMatrix3d fineAxes = cRefineOBBAxes(pts, perturbed, fine, &fineStats);
        double fineErr = cMax(axisAngleDeg(fineAxes.getCol0(), R.getCol0()),
                         cMax(axisAngleDeg(fineAxes.getCol1(), R.getCol1()),
                              axisAngleDeg(fineAxes.getCol2(), R.getCol2())));
        CHECK(fineStats.m_finalVolume <= stats.m_finalVolume);
        CHECK(fineErr < 0.2);
        printf("perturbed frame, steps {5,1,0.2} deg: volume %.4f -> %.4f (exact 8), axis error %.3f deg, "
               "%d iterations, %d rotations, %d evaluations\n",
               fineStats.m_initialVolume, fineStats.m_finalVolume, fineErr,
               fineStats.m_iterations, fineStats.m_acceptedRotations, fineStats.m_evaluations);

        // the refinement keeps an orthonormal right-handed frame
        CHECK(approx(refined * cTranspose(refined), cIdentity3d(), 1e-12));
        CHECK(approx(refined.det(), 1.0, 1e-12));
    }

    //--------------------------------------------------------------------------
    // an optimal PCA frame is left unchanged (no rotation accepted on noise)
    //--------------------------------------------------------------------------
    {
        const cMatrix3d R = rotationDeg(cVector3d(-1, 0.5, 2), 71);
        std::vector<cVector3d> pts = boxCorners(cVector3d(2, 1, 0.5), R, cVector3d(0, 0, 0));
        cOBBRefineStats stats;
        cCollisionOBBBox box = cComputeOBB(pts, NULL, cOBBRefineSettings(), &stats);
        CHECK(stats.m_acceptedRotations == 0);
        CHECK(stats.m_iterations == 1);
        CHECK(approx(box.getHalfExtents(), cVector3d(2, 1, 0.5), 1e-5));

        // disabled refinement reports the PCA volume
        cOBBRefineStats off;
        cComputeOBB(pts, NULL, noRefinement(), &off);
        CHECK(off.m_iterations == 0 && off.m_acceptedRotations == 0);
        CHECK(approx(off.m_initialVolume, off.m_finalVolume, 0.0 + 1e-300));
    }

    //--------------------------------------------------------------------------
    // random point clouds: the refined box is never larger than the PCA box
    //--------------------------------------------------------------------------
    {
        std::normal_distribution<double> gauss(0.0, 1.0);
        double sumGain = 0.0;
        double maxGain = 0.0;
        int improved = 0;
        const int trials = 30;
        for (int trial=0; trial<trials; trial++)
        {
            // anisotropic, sheared gaussian cloud
            cMatrix3d shape = randomRotation(rng) *
                              cMatrix3d(3, 0.8, 0,  0, 1.5, 0.4,  0, 0, 0.7) *
                              randomRotation(rng);
            std::vector<cVector3d> pts;
            for (int i=0; i<60; i++)
                pts.push_back(shape * randomVector(rng, gauss));

            cOBBRefineStats pca, ref;
            cCollisionOBBBox boxPCA = cComputeOBB(pts, NULL, noRefinement(), &pca);
            cCollisionOBBBox boxRef = cComputeOBB(pts, NULL, cOBBRefineSettings(), &ref);

            CHECK(ref.m_finalVolume <= pca.m_finalVolume * (1.0 + 1e-12));
            CHECK(containsAll(boxPCA, pts));
            CHECK(containsAll(boxRef, pts));

            double gain = 1.0 - ref.m_finalVolume / pca.m_finalVolume;
            sumGain += gain;
            maxGain = cMax(maxGain, gain);
            if (gain > 1e-9) improved++;
        }
        printf("random clouds: refinement improved %d/%d boxes, mean volume reduction %.2f%%, max %.2f%%\n",
               improved, trials, 100.0 * sumGain / trials, 100.0 * maxGain);
        CHECK(improved > trials / 2);
    }

    //--------------------------------------------------------------------------
    // degenerate eigenvalues: rotated cube, PCA axes are arbitrary
    //--------------------------------------------------------------------------
    {
        const cMatrix3d R = randomRotation(rng);
        std::vector<cVector3d> pts = boxCorners(cVector3d(1, 1, 1), R, cVector3d(0, 0, 0));

        cOBBRefineStats pca, ref, wide;
        cComputeOBB(pts, NULL, noRefinement(), &pca);
        cComputeOBB(pts, NULL, cOBBRefineSettings(), &ref);

        cOBBRefineSettings wideSettings;
        wideSettings.m_maxIterations = 30;
        wideSettings.m_anglesDeg.clear();
        wideSettings.m_anglesDeg.push_back(15.0);
        wideSettings.m_anglesDeg.push_back(5.0);
        wideSettings.m_anglesDeg.push_back(1.0);
        wideSettings.m_anglesDeg.push_back(0.2);
        cComputeOBB(pts, NULL, wideSettings, &wide);

        CHECK(ref.m_finalVolume <= pca.m_finalVolume * (1.0 + 1e-12));
        CHECK(wide.m_finalVolume <= pca.m_finalVolume * (1.0 + 1e-12));
        printf("cube (exact volume 8): PCA %.4f, refined 5 x {5,1} deg %.4f, refined 30 x {15,5,1,0.2} deg %.4f\n",
               pca.m_finalVolume, ref.m_finalVolume, wide.m_finalVolume);
    }

    //--------------------------------------------------------------------------
    // flat acute triangle: the minimum bounding rectangle has area 2 x 4 = 8
    // and is aligned with an edge; the refinement minimizes the area
    //--------------------------------------------------------------------------
    {
        const cMatrix3d R = rotationDeg(cVector3d(1, 1, 0), 40);
        std::vector<cVector3d> pts;
        pts.push_back(R * cVector3d(0, 0, 0));
        pts.push_back(R * cVector3d(4, 0, 0));
        pts.push_back(R * cVector3d(1, 2, 0));

        cCollisionOBBBox boxPCA = cComputeOBB(pts, NULL, noRefinement());
        cCollisionOBBBox boxRef = cComputeOBB(pts, NULL);
        cVector3d hP = boxPCA.getHalfExtents();
        cVector3d hR = boxRef.getHalfExtents();
        double areaPCA = 4.0 * hP(0) * hP(1);
        double areaRef = 4.0 * hR(0) * hR(1);

        CHECK(hR(2) < 1e-6);
        CHECK(areaRef < areaPCA);
        CHECK(areaRef < 8.0 * 1.05);
        CHECK(containsAll(boxRef, pts));
        printf("flat triangle (min rectangle area 8): PCA %.4f, refined %.4f\n", areaPCA, areaRef);
    }

    //--------------------------------------------------------------------------
    // cost: refinement of a hull with many vertices
    //--------------------------------------------------------------------------
    {
        std::uniform_real_distribution<double> uni(-1.0, 1.0);
        std::vector<cVector3d> pts;
        while (pts.size() < 20000)
        {
            cVector3d p = randomVector(rng, uni);
            if (p.length() > 0.1) pts.push_back(cMatrix3d(2, 0, 0, 0, 1, 0, 0, 0, 0.5) * cNormalize(p));
        }

        cConvexHull hull;
        cPrecisionClock clock;
        clock.start(true);
        cOBBRefineStats stats;
        cComputeOBB(pts, &hull, cOBBRefineSettings(), &stats);
        double total = clock.getCurrentTimeSeconds();
        printf("ellipsoid, %u hull vertices: full fit %.3f s, %d evaluations\n",
               hull.getNumVertices(), total, stats.m_evaluations);
    }

    return (testResult());
}
