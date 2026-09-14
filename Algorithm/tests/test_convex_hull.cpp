//==============================================================================
/*
    Tests for cConvexHull (Quickhull).
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CConvexHull.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
#include "math/CMaths.h"
#include "timers/CPrecisionClock.h"
//------------------------------------------------------------------------------
#include <map>
#include <random>
#include <utility>
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Checks that a 3D hull is a closed, consistently oriented, convex polyhedron
// that contains every input point.
//------------------------------------------------------------------------------
static bool isValidHull(const cConvexHull& a_hull, const std::vector<cVector3d>& a_points, double a_eps)
{
    const std::vector<cVector3d>& V = a_hull.getVertices();
    const std::vector<unsigned int>& T = a_hull.getTriangles();
    if (a_hull.getDimension() != 3 || T.empty()) return (false);

    // every directed edge appears once, and its reverse appears once
    std::map<std::pair<unsigned int, unsigned int>, int> edges;
    for (size_t t=0; t<T.size(); t+=3)
    {
        for (int k=0; k<3; k++)
        {
            edges[std::make_pair(T[t+k], T[t+(k+1)%3])]++;
        }
    }
    for (auto it=edges.begin(); it!=edges.end(); ++it)
    {
        if (it->second != 1) return (false);
        if (edges.count(std::make_pair(it->first.second, it->first.first)) != 1) return (false);
    }

    // Euler characteristic of a sphere: V - E + F = 2
    long numV = (long)V.size();
    long numE = (long)edges.size() / 2;
    long numF = (long)T.size() / 3;
    if (numV - numE + numF != 2) return (false);

    // convexity: every input point lies behind every face
    for (size_t t=0; t<T.size(); t+=3)
    {
        const cVector3d& a = V[T[t]];
        const cVector3d& b = V[T[t+1]];
        const cVector3d& c = V[T[t+2]];
        cVector3d n = cCross(b - a, c - a);
        if (n.length() <= 0.0) return (false);
        n.normalize();
        for (size_t i=0; i<a_points.size(); i++)
        {
            if (cDot(n, a_points[i] - a) > a_eps) return (false);
        }
    }

    return (true);
}

//------------------------------------------------------------------------------
// Points on the surface of an axis-aligned box, each face tessellated as a
// regular grid (many coplanar vertices).
//------------------------------------------------------------------------------
static void addTessellatedBox(std::vector<cVector3d>& a_points, const cVector3d& a_half, int a_grid)
{
    for (int axis=0; axis<3; axis++)
    {
        int u = (axis + 1) % 3;
        int v = (axis + 2) % 3;
        for (int side=-1; side<=1; side+=2)
        {
            for (int i=0; i<=a_grid; i++)
            {
                for (int j=0; j<=a_grid; j++)
                {
                    cVector3d p;
                    p(axis) = side * a_half(axis);
                    p(u) = a_half(u) * (-1.0 + 2.0 * i / a_grid);
                    p(v) = a_half(v) * (-1.0 + 2.0 * j / a_grid);
                    a_points.push_back(p);
                }
            }
        }
    }
}


int main()
{
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    cConvexHull hull;

    //--------------------------------------------------------------------------
    // degenerate inputs
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> pts;
        CHECK(!hull.compute(pts));
        CHECK(hull.getDimension() == -1);

        // coincident points
        pts.assign(5, cVector3d(1, 2, 3));
        CHECK(hull.compute(pts));
        CHECK(hull.getDimension() == 0);
        CHECK(hull.getNumVertices() == 1);

        // collinear points: the two extremes are kept
        pts.clear();
        cVector3d dir = cNormalize(cVector3d(1, -2, 0.5));
        for (int i=0; i<50; i++) pts.push_back(cVector3d(1, 1, 1) + uni(rng) * dir);
        pts.push_back(cVector3d(1, 1, 1) + 3.0 * dir);
        pts.push_back(cVector3d(1, 1, 1) - 2.0 * dir);
        CHECK(hull.compute(pts));
        CHECK(hull.getDimension() == 1);
        CHECK(hull.getNumVertices() == 2 && hull.getNumTriangles() == 0);
        CHECK(approx(cDistance(hull.getVertices()[0], hull.getVertices()[1]), 5.0, 1e-12));

        // coplanar points: tessellated square in a tilted plane
        pts.clear();
        cMatrix3d rot;
        rot.setAxisAngleRotationDeg(cVector3d(1, 1, 0), 30);
        for (int i=0; i<=20; i++)
            for (int j=0; j<=20; j++)
                pts.push_back(rot * cVector3d(i / 10.0 - 1.0, j / 10.0 - 1.0, 0.0));
        CHECK(hull.compute(pts));
        CHECK(hull.getDimension() == 2);
        CHECK(hull.getNumVertices() == 4);
        CHECK(hull.getNumTriangles() == 2);
    }

    //--------------------------------------------------------------------------
    // tessellated box with interior points: only the 8 corners remain
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> pts;
        cVector3d half(2.0, 1.0, 0.5);
        addTessellatedBox(pts, half, 20);
        for (int i=0; i<1000; i++)
            pts.push_back(randomInBox(rng, half));

        CHECK(hull.compute(pts));
        CHECK(hull.getNumVertices() == 8);
        CHECK(hull.getNumTriangles() == 12);
        CHECK(isValidHull(hull, pts, 1e-12));
        printf("tessellated box: %d points -> %u hull vertices, %u triangles\n",
               (int)pts.size(), hull.getNumVertices(), hull.getNumTriangles());

        // same box far from the origin, with duplicated points
        // note: "far" cannot be used as a name on Windows (empty macro in minwindef.h)
        std::vector<cVector3d> shifted;
        for (size_t i=0; i<pts.size(); i++)
        {
            shifted.push_back(pts[i] + cVector3d(1e4, -2e4, 3e4));
            if (i % 7 == 0) shifted.push_back(shifted.back());
        }
        CHECK(hull.compute(shifted));
        CHECK(hull.getNumVertices() == 8);
        CHECK(isValidHull(hull, shifted, 1e-8));
    }

    //--------------------------------------------------------------------------
    // random points in a ball, and points that all lie on a sphere
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> pts;
        while (pts.size() < 5000)
        {
            cVector3d p = randomVector(rng, uni);
            if (p.length() <= 1.0) pts.push_back(p);
        }
        CHECK(hull.compute(pts));
        CHECK(isValidHull(hull, pts, 1e-12));
        printf("ball: %d points -> %u hull vertices\n", (int)pts.size(), hull.getNumVertices());

        pts.clear();
        while (pts.size() < 2000)
        {
            cVector3d p = randomVector(rng, uni);
            if (p.length() > 0.1) pts.push_back(cNormalize(p));
        }
        CHECK(hull.compute(pts));
        CHECK(isValidHull(hull, pts, 1e-12));
        CHECK(hull.getNumVertices() == 2000);
    }

    //--------------------------------------------------------------------------
    // performance: 100000 points on a sphere (every point is a hull vertex)
    //--------------------------------------------------------------------------
    {
        std::vector<cVector3d> pts;
        while (pts.size() < 100000)
        {
            cVector3d p = randomVector(rng, uni);
            if (p.length() > 0.1) pts.push_back(cNormalize(p));
        }
        cPrecisionClock clock;
        clock.start(true);
        CHECK(hull.compute(pts));
        double elapsed = clock.getCurrentTimeSeconds();
        printf("sphere: %d points -> %u hull vertices in %.3f s\n",
               (int)pts.size(), hull.getNumVertices(), elapsed);
        CHECK(hull.getNumVertices() > 99000);
    }

    return (testResult());
}
