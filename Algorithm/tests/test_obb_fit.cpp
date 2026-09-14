//==============================================================================
/*
    Tests for the tight OBB fit: convex hull preprocessing + continuous PCA.
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBBFit.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
#include <algorithm>
#include <random>
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Appends triangle (a, b, c) subdivided into a_n * a_n equal triangles.
//------------------------------------------------------------------------------
static void addSubdividedTriangle(std::vector<cVector3d>& a_vertices,
                                  std::vector<unsigned int>& a_triangles,
                                  const cVector3d& a, const cVector3d& b, const cVector3d& c,
                                  int a_n)
{
    std::vector<unsigned int> index((a_n + 1) * (a_n + 1), 0);
    for (int i=0; i<=a_n; i++)
    {
        for (int j=0; i+j<=a_n; j++)
        {
            index[i * (a_n + 1) + j] = (unsigned int)a_vertices.size();
            a_vertices.push_back(a + (double)i / a_n * (b - a) + (double)j / a_n * (c - a));
        }
    }
    auto id = [&](int i, int j) { return index[i * (a_n + 1) + j]; };
    for (int i=0; i<a_n; i++)
    {
        for (int j=0; i+j<a_n; j++)
        {
            a_triangles.push_back(id(i, j));
            a_triangles.push_back(id(i + 1, j));
            a_triangles.push_back(id(i, j + 1));
            if (i + j + 1 < a_n)
            {
                a_triangles.push_back(id(i + 1, j));
                a_triangles.push_back(id(i + 1, j + 1));
                a_triangles.push_back(id(i, j + 1));
            }
        }
    }
}

//------------------------------------------------------------------------------
// Angle (degrees) between a direction and an axis, ignoring the sign.
//------------------------------------------------------------------------------
static double axisAngleDeg(const cVector3d& a_dir, const cVector3d& a_axis)
{
    double c = cClamp(fabs(cDot(cNormalize(a_dir), cNormalize(a_axis))), 0.0, 1.0);
    return (cRadToDeg(acos(c)));
}

static double volume(const cCollisionOBBBox& a_box)
{
    cVector3d h = a_box.getHalfExtents();
    return (8.0 * h(0) * h(1) * h(2));
}

static bool containsAll(const cCollisionOBBBox& a_box, const std::vector<cVector3d>& a_points)
{
    for (size_t i=0; i<a_points.size(); i++)
    {
        if (!a_box.contains(a_points[i])) return (false);
    }
    return (true);
}


int main()
{
    std::mt19937 rng(2024);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    //--------------------------------------------------------------------------
    // continuous covariance of the unit square: mean (0.5, 0.5), var = 1/12
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> V;
        V.push_back(cVector3d(0, 0, 0));
        V.push_back(cVector3d(1, 0, 0));
        V.push_back(cVector3d(1, 1, 0));
        V.push_back(cVector3d(0, 1, 0));
        unsigned int tri[] = { 0, 1, 2,  0, 2, 3 };
        std::vector<unsigned int> T(tri, tri + 6);

        cVector3d mu;
        cMatrix3d C;
        double area = 0.0;
        CHECK(cComputeSurfaceCovariance(V, T, mu, C, &area));
        CHECK(approx(area, 1.0, 1e-14));
        CHECK(approx(mu, cVector3d(0.5, 0.5, 0.0), 1e-14));
        CHECK(approx(C, cMatrix3d(1.0/12, 0, 0,  0, 1.0/12, 0,  0, 0, 0), 1e-14));
    }

    //--------------------------------------------------------------------------
    // tessellation bias: a 4 x 1 rectangle made of two triangles, one of them
    // subdivided into 4096 small triangles
    //--------------------------------------------------------------------------
    {
        const cVector3d p00(0, 0, 0), p40(4, 0, 0), p41(4, 1, 0), p01(0, 1, 0);

        std::vector<cVector3d> coarseV, fineV;
        std::vector<unsigned int> coarseT, fineT;
        addSubdividedTriangle(coarseV, coarseT, p00, p40, p41, 1);
        addSubdividedTriangle(coarseV, coarseT, p00, p41, p01, 1);
        addSubdividedTriangle(fineV, fineT, p00, p40, p41, 1);
        addSubdividedTriangle(fineV, fineT, p00, p41, p01, 64);

        cVector3d muCoarse, muFine;
        cMatrix3d cCoarse, cFine;
        CHECK(cComputeSurfaceCovariance(coarseV, coarseT, muCoarse, cCoarse));
        CHECK(cComputeSurfaceCovariance(fineV, fineT, muFine, cFine));

        // continuous PCA does not depend on the tessellation
        CHECK(approx(muFine, cVector3d(2.0, 0.5, 0.0), 1e-12));
        CHECK(approx(muFine, muCoarse, 1e-12));
        CHECK(approx(cFine, cCoarse, 1e-12));
        cMatrix3d axesFine = cComputePrincipalAxes(cFine);
        double angleContinuous = axisAngleDeg(axesFine.getCol0(), cVector3d(1, 0, 0));
        CHECK(angleContinuous < 1e-9);

        // vertex-based PCA is pulled toward the dense triangle
        cVector3d muVertex;
        cMatrix3d cVertex;
        cComputePointCovariance(fineV, muVertex, cVertex);
        double angleVertex = axisAngleDeg(cComputePrincipalAxes(cVertex).getCol0(), cVector3d(1, 0, 0));
        CHECK(cDistance(muVertex, muFine) > 0.5);
        CHECK(angleVertex > 5.0);

        printf("tessellation bias (%d triangles): centroid error vertex PCA = %.3f, continuous = %.1e\n",
               (int)fineT.size() / 3, cDistance(muVertex, muFine), cDistance(muFine, cVector3d(2.0, 0.5, 0.0)));
        printf("                                  major axis tilt  vertex PCA = %.2f deg, continuous = %.1e deg\n",
               angleVertex, angleContinuous);
    }

    //--------------------------------------------------------------------------
    // rotated box sampled with an asymmetric dense patch and interior points
    //--------------------------------------------------------------------------
    {
        const cVector3d half(2.0, 1.0, 0.5);
        cMatrix3d R;
        R.setAxisAngleRotationDeg(cVector3d(1, 2, 3), 37);
        const cVector3d t(10, -5, 3);

        std::vector<cVector3d> local;
        for (int i=0; i<8; i++)
            local.push_back(cVector3d((i & 1) ? half(0) : -half(0),
                                      (i & 2) ? half(1) : -half(1),
                                      (i & 4) ? half(2) : -half(2)));
        for (int i=0; i<=40; i++)                       // dense patch on one corner of the +X face
            for (int j=0; j<=40; j++)
                local.push_back(cVector3d(half(0), i / 40.0 * half(1), j / 40.0 * half(2)));
        for (int i=0; i<2000; i++)
            local.push_back(randomInBox(rng, half));

        std::vector<cVector3d> pts;
        for (size_t i=0; i<local.size(); i++) pts.push_back(R * local[i] + t);

        cConvexHull hull;
        cCollisionOBBBox box = cComputeOBB(pts, &hull);
        CHECK(hull.getNumVertices() == 8);
        CHECK(approx(box.getHalfExtents(), half, 1e-5));
        CHECK(approx(box.getCenter(), t, 1e-5));
        CHECK(axisAngleDeg(box.getAxis(0), R.getCol0()) < 1e-3);
        CHECK(axisAngleDeg(box.getAxis(1), R.getCol1()) < 1e-3);
        CHECK(axisAngleDeg(box.getAxis(2), R.getCol2()) < 1e-3);
        CHECK(containsAll(box, pts));

        // vertex-based PCA on the same points gives a looser box
        cVector3d mean;
        cMatrix3d cov;
        cComputePointCovariance(pts, mean, cov);
        cCollisionOBBBox vertexBox = cComputeOBBFromAxes(pts, cComputePrincipalAxes(cov));
        CHECK(containsAll(vertexBox, pts));
        CHECK(volume(box) < volume(vertexBox));
        printf("rotated box (%d points): volume exact = 8, continuous = %.5f, vertex PCA = %.5f\n",
               (int)pts.size(), volume(box), volume(vertexBox));
    }

    //--------------------------------------------------------------------------
    // random poses, including far from the origin (float storage)
    //--------------------------------------------------------------------------
    {
        std::uniform_real_distribution<double> size(0.1, 3.0);
        int failsBefore = g_fails;
        for (int trial=0; trial<50; trial++)
        {
            double e[3] = { size(rng), size(rng), size(rng) };
            std::sort(e, e + 3);
            if ((e[1] - e[0] < 0.05) || (e[2] - e[1] < 0.05)) { trial--; continue; }
            const cVector3d half(e[2], e[1], e[0]);
            const cMatrix3d R = randomRotation(rng);
            const double scale = (trial < 25) ? 10.0 : 1e5;
            const cVector3d t = scale * randomVector(rng, uni);

            std::vector<cVector3d> pts;
            for (int i=0; i<8; i++)
                pts.push_back(R * cVector3d((i & 1) ? half(0) : -half(0),
                                            (i & 2) ? half(1) : -half(1),
                                            (i & 4) ? half(2) : -half(2)) + t);
            for (int i=0; i<200; i++)
                pts.push_back(R * randomInBox(rng, half) + t);

            cCollisionOBBBox box = cComputeOBB(pts);
            // at 1e5 the float center can be off by half a float spacing (0.0039)
            // per coordinate, which the half-extents must absorb
            double tolerance = (trial < 25) ? 1e-5 : 1.2e-2;
            CHECK(containsAll(box, pts));
            CHECK(approx(box.getHalfExtents(), half, tolerance));
            if (!approx(box.getHalfExtents(), half, tolerance))
            {
                printf("  trial %d: expected %s, got %s\n", trial,
                       half.str(6).c_str(), box.getHalfExtents().str(6).c_str());
            }
        }
        printf("random poses: %s\n", (g_fails == failsBefore) ? "50/50 tight and enclosing" : "FAILED");
    }

    //--------------------------------------------------------------------------
    // CHAI3D mesh: box primitive with position and rotation
    //--------------------------------------------------------------------------
    {
        cMatrix3d rot;
        rot.setAxisAngleRotationDeg(cVector3d(0.3, -1, 0.5), 50);
        const cVector3d pos(1, 2, 3);

        cMesh* mesh = new cMesh();
        cCreateBox(mesh, 4.0, 2.0, 1.0, pos, rot);

        cConvexHull hull;
        cCollisionOBBBox box = cComputeOBB(mesh, &hull);
        CHECK(hull.getNumVertices() == 8);
        CHECK(approx(box.getHalfExtents(), cVector3d(2.0, 1.0, 0.5), 1e-5));
        CHECK(approx(box.getCenter(), pos, 1e-5));
        CHECK(axisAngleDeg(box.getAxis(0), rot.getCol0()) < 1e-3);
        CHECK(axisAngleDeg(box.getAxis(2), rot.getCol2()) < 1e-3);

        std::vector<cVector3d> verts;
        for (unsigned int i=0; i<mesh->getNumVertices(); i++) verts.push_back(mesh->m_vertices->getLocalPos(i));
        CHECK(containsAll(box, verts));
        printf("cMesh box: %u vertices, %u triangles -> OBB half-extents %s\n",
               mesh->getNumVertices(), mesh->getNumTriangles(), box.getHalfExtents().str(4).c_str());
        delete mesh;
    }

    //--------------------------------------------------------------------------
    // degenerate inputs
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> pts;
        cCollisionOBBBox box = cComputeOBB(pts);
        CHECK(approx(box.getHalfExtents(), cVector3d(0, 0, 0), 0.0 + 1e-30));

        pts.assign(3, cVector3d(5, 5, 5));
        box = cComputeOBB(pts);
        CHECK(approx(box.getCenter(), cVector3d(5, 5, 5), 1e-6));
        CHECK(approx(box.getHalfExtents(), cVector3d(0, 0, 0), 1e-6));
        CHECK(containsAll(box, pts));

        // segment from (0,0,0) to (2,2,0)
        pts.clear();
        for (int i=0; i<=10; i++) pts.push_back(cVector3d(0.2 * i, 0.2 * i, 0.0));
        box = cComputeOBB(pts);
        CHECK(approx(box.getHalfExtents(), cVector3d(sqrt(2.0), 0, 0), 1e-5));
        CHECK(axisAngleDeg(box.getAxis(0), cVector3d(1, 1, 0)) < 1e-3);
        CHECK(containsAll(box, pts));

        // flat 3 x 1 rectangle in a tilted plane
        pts.clear();
        cMatrix3d R;
        R.setAxisAngleRotationDeg(cVector3d(1, 0, 1), 25);
        for (int i=0; i<=30; i++)
            for (int j=0; j<=10; j++)
                pts.push_back(R * cVector3d(-1.5 + 0.1 * i, -0.5 + 0.1 * j, 0.0));
        box = cComputeOBB(pts);
        CHECK(approx(box.getHalfExtents(), cVector3d(1.5, 0.5, 0.0), 1e-5));
        CHECK(axisAngleDeg(box.getAxis(2), R.getCol2()) < 1e-3);
        CHECK(containsAll(box, pts));
    }

    return (testResult());
}
