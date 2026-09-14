//==============================================================================
/*
    Tests for the intersection primitives of the OBB tree. Each test is
    compared with an independent formulation on random configurations.
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBBIntersect.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
#include <cfloat>
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Reference triangle-triangle test for triangles in general position: they
// intersect if and only if an edge of one of them crosses the other one.
//------------------------------------------------------------------------------
static bool triTriByEdges(const cVector3d P[3], const cVector3d Q[3])
{
    double t;
    for (int i=0; i<3; i++)
    {
        if (cTestSegmentTriangle(P[i], P[(i+1)%3], Q[0], Q[1], Q[2], t)) return (true);
        if (cTestSegmentTriangle(Q[i], Q[(i+1)%3], P[0], P[1], P[2], t)) return (true);
    }
    return (false);
}

//------------------------------------------------------------------------------
// Corners and the 12 face triangles of a box.
//------------------------------------------------------------------------------
static void boxGeometry(const cOBBShape& a_box, cVector3d a_corners[8], int a_faces[12][3])
{
    for (int i=0; i<8; i++)
    {
        cVector3d local((i & 1) ? a_box.m_half(0) : -a_box.m_half(0),
                        (i & 2) ? a_box.m_half(1) : -a_box.m_half(1),
                        (i & 4) ? a_box.m_half(2) : -a_box.m_half(2));
        a_corners[i] = a_box.m_rot * local + a_box.m_center;
    }
    static const int faces[12][3] = { {0,1,3},{0,3,2}, {4,6,7},{4,7,5}, {0,4,5},{0,5,1},
                                      {2,3,7},{2,7,6}, {0,2,6},{0,6,4}, {1,5,7},{1,7,3} };
    for (int f=0; f<12; f++)
        for (int k=0; k<3; k++)
            a_faces[f][k] = faces[f][k];
}

static bool insideBox(const cOBBShape& a_box, const cVector3d& a_point)
{
    const cVector3d d = a_point - a_box.m_center;
    return ((fabs(cDot(d, a_box.m_rot.getCol0())) <= a_box.m_half(0)) &&
            (fabs(cDot(d, a_box.m_rot.getCol1())) <= a_box.m_half(1)) &&
            (fabs(cDot(d, a_box.m_rot.getCol2())) <= a_box.m_half(2)));
}

//------------------------------------------------------------------------------
// Reference box-box test: faces intersect, or one box contains the other.
//------------------------------------------------------------------------------
static bool boxBoxByFaces(const cOBBShape& a_a, const cOBBShape& a_b)
{
    cVector3d ca[8], cb[8];
    int fa[12][3], fb[12][3];
    boxGeometry(a_a, ca, fa);
    boxGeometry(a_b, cb, fb);
    if (insideBox(a_b, ca[0]) || insideBox(a_a, cb[0])) return (true);
    for (int i=0; i<12; i++)
        for (int j=0; j<12; j++)
            if (cTestTriangleTriangle(ca[fa[i][0]], ca[fa[i][1]], ca[fa[i][2]],
                                      cb[fb[j][0]], cb[fb[j][1]], cb[fb[j][2]])) return (true);
    return (false);
}

static cOBBShape randomBox(std::mt19937& a_rng, double a_spread)
{
    std::uniform_real_distribution<double> pos(-a_spread, a_spread);
    std::uniform_real_distribution<double> size(0.05, 1.0);
    cOBBShape box;
    box.m_center = randomVector(a_rng, pos);
    box.m_rot = randomRotation(a_rng);
    box.m_half = randomVector(a_rng, size);
    return (box);
}


int main()
{
    std::mt19937 rng(99);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    //--------------------------------------------------------------------------
    // triangle-triangle: hand-made cases
    //--------------------------------------------------------------------------
    {
        const cVector3d a(0, 0, 0), b(1, 0, 0), c(0, 1, 0);

        // crossing
        CHECK( cTestTriangleTriangle(a, b, c, cVector3d(0.25, 0.25, -1), cVector3d(0.25, 0.25, 1), cVector3d(2, 2, 0)));
        // parallel planes
        CHECK(!cTestTriangleTriangle(a, b, c, cVector3d(0, 0, 0.5), cVector3d(1, 0, 0.5), cVector3d(0, 1, 0.5)));
        // shared vertex
        CHECK( cTestTriangleTriangle(a, b, c, cVector3d(1, 0, 0), cVector3d(2, 0, 1), cVector3d(2, 1, -1)));
        // coplanar, overlapping
        CHECK( cTestTriangleTriangle(a, b, c, cVector3d(0.2, 0.2, 0), cVector3d(1.2, 0.2, 0), cVector3d(0.2, 1.2, 0)));
        // coplanar, disjoint (only an in-plane edge normal separates them)
        CHECK(!cTestTriangleTriangle(a, b, c, cVector3d(1, 1, 0), cVector3d(2, 1, 0), cVector3d(1, 2, 0)));
        // crosses the plane of the first triangle outside of it
        CHECK(!cTestTriangleTriangle(a, b, c, cVector3d(0.8, 0.8, -1), cVector3d(0.8, 0.8, 1), cVector3d(2, 0.8, 0)));
    }

    //--------------------------------------------------------------------------
    // triangle-triangle: random triangles against the edge formulation
    //--------------------------------------------------------------------------
    {
        int mismatches = 0, hits = 0;
        const int trials = 20000;
        for (int i=0; i<trials; i++)
        {
            cVector3d P[3], Q[3];
            for (int k=0; k<3; k++) P[k] = randomVector(rng, uni);
            const cVector3d offset = 0.8 * randomVector(rng, uni);
            for (int k=0; k<3; k++) Q[k] = randomVector(rng, uni) + offset;

            const bool sat = cTestTriangleTriangle(P[0], P[1], P[2], Q[0], Q[1], Q[2]);
            const bool ref = triTriByEdges(P, Q);
            if (sat != ref) mismatches++;
            if (ref) hits++;
        }
        CHECK(mismatches == 0);
        printf("triangle-triangle: %d random pairs (%d intersecting), %d mismatches\n", trials, hits, mismatches);
    }

    //--------------------------------------------------------------------------
    // box-box: random boxes against the face formulation
    //--------------------------------------------------------------------------
    {
        int mismatches = 0, hits = 0;
        const int trials = 5000;
        for (int i=0; i<trials; i++)
        {
            const cOBBShape A = randomBox(rng, 1.5);
            const cOBBShape B = randomBox(rng, 1.5);
            const bool sat = cTestOBBOBB(A, B);
            const bool ref = boxBoxByFaces(A, B);
            if (sat != ref) mismatches++;
            if (ref) hits++;
        }
        CHECK(mismatches == 0);
        printf("box-box:           %d random pairs (%d overlapping), %d mismatches\n", trials, hits, mismatches);
    }

    //--------------------------------------------------------------------------
    // segment-box: random segments against the face formulation
    //--------------------------------------------------------------------------
    {
        int mismatches = 0, hits = 0, badEntry = 0;
        const int trials = 20000;
        for (int i=0; i<trials; i++)
        {
            const cOBBShape box = randomBox(rng, 0.5);
            const cVector3d a = 2.0 * randomVector(rng, uni);
            const cVector3d b = 2.0 * randomVector(rng, uni);

            cVector3d corners[8];
            int faces[12][3];
            boxGeometry(box, corners, faces);
            double tFaces = DBL_MAX, t;
            for (int f=0; f<12; f++)
                if (cTestSegmentTriangle(a, b, corners[faces[f][0]], corners[faces[f][1]], corners[faces[f][2]], t))
                    tFaces = cMin(tFaces, t);
            const bool aInside = insideBox(box, a);
            const bool ref = aInside || insideBox(box, b) || (tFaces < DBL_MAX);

            double tEnter = -1.0;
            const bool slab = cTestSegmentOBB(a, b, box, 0.0, tEnter);
            if (slab != ref) mismatches++;
            if (ref) hits++;
            if (slab && ref && !approx(tEnter, aInside ? 0.0 : tFaces, 1e-9)) badEntry++;
        }
        CHECK(mismatches == 0);
        CHECK(badEntry == 0);
        printf("segment-box:       %d random segments (%d hitting), %d mismatches, %d wrong entry points\n",
               trials, hits, mismatches, badEntry);
    }

    //--------------------------------------------------------------------------
    // segment-triangle: against CHAI3D's cIntersectionSegmentTriangle
    //--------------------------------------------------------------------------
    {
        int mismatches = 0, hits = 0, badPoint = 0;
        const int trials = 20000;
        for (int i=0; i<trials; i++)
        {
            cVector3d v[3];
            for (int k=0; k<3; k++) v[k] = randomVector(rng, uni);
            const cVector3d a = 1.5 * randomVector(rng, uni);
            const cVector3d b = 1.5 * randomVector(rng, uni);

            double t;
            const bool ours = cTestSegmentTriangle(a, b, v[0], v[1], v[2], t);
            cVector3d point, normal;
            double p01, p02;
            const bool chai = cIntersectionSegmentTriangle(a, b, v[0], v[1], v[2], true, true, point, normal, p01, p02);
            if (ours != chai) mismatches++;
            if (chai) hits++;
            if (ours && chai && !approx(a + t * (b - a), point, 1e-9)) badPoint++;
        }
        CHECK(mismatches == 0);
        CHECK(badPoint == 0);
        printf("segment-triangle:  %d random segments (%d hitting), %d mismatches with CHAI3D, %d different points\n",
               trials, hits, mismatches, badPoint);
    }

    return (testResult());
}
