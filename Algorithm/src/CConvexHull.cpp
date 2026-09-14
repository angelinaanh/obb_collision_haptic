//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CConvexHull.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CConvexHull.h"
//------------------------------------------------------------------------------
#include "math/CMaths.h"
//------------------------------------------------------------------------------
#include <algorithm>
#include <cfloat>
#include <cmath>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace {
//------------------------------------------------------------------------------

//! Triangular face of the hull under construction.
struct cHullFace
{
    //! Point indices, counter clockwise seen from outside.
    int m_v[3];

    //! Neighbor face across edge (m_v[k], m_v[(k+1)%3]).
    int m_nb[3];

    //! Unit outward normal.
    cVector3d m_normal;

    //! Plane offset: dot(m_normal, x) = m_offset for points x on the plane.
    double m_offset;

    //! Indices of the points located in front of this face.
    std::vector<int> m_outside;

    //! __false__ once the face has been replaced.
    bool m_alive;

    //! Visibility of the face from the current eye point.
    bool m_visible;

    //! Stamp of the last visibility search that visited this face.
    unsigned int m_visit;

    //! Signed distance from the plane of the face to a point.
    inline double distance(const cVector3d& a_point) const
    {
        return (cDot(m_normal, a_point) - m_offset);
    }
};


//! Edge between a visible face and a non-visible face.
struct cHorizonEdge
{
    int m_a;
    int m_b;
    int m_face;     // non-visible face across the edge
};


//! Selection of the farthest point with a tie-break.
/*!
    A candidate replaces the current best one if it is farther by more than
    the tolerance, or as far (within the tolerance) and larger on the
    secondary criterion. Points that tie on distance typically form a flat
    patch, e.g. all the vertices of a tessellated face parallel to the
    reference plane. Choosing the one farthest from a fixed point selects an
    extreme point of that patch, i.e. a true hull vertex, instead of a vertex
    lying inside a flat face or on a straight edge.
*/
struct cFarthest
{
    int m_index;
    double m_dist;
    double m_second;

    cFarthest() : m_index(-1), m_dist(-DBL_MAX), m_second(-DBL_MAX) {}

    void test(int a_index, double a_dist, double a_second, double a_tolerance)
    {
        if (a_dist > m_dist + a_tolerance)
        {
            m_index = a_index;
            m_dist = a_dist;
            m_second = a_second;
        }
        else if ((a_dist >= m_dist - a_tolerance) && (a_second > m_second))
        {
            m_index = a_index;
            m_dist = cMax(m_dist, a_dist);
            m_second = a_second;
        }
    }
};


//! Creates a face and computes its plane.
cHullFace newFace(const std::vector<cVector3d>& a_points, int a_v0, int a_v1, int a_v2)
{
    cHullFace face;
    face.m_v[0] = a_v0;
    face.m_v[1] = a_v1;
    face.m_v[2] = a_v2;
    face.m_nb[0] = face.m_nb[1] = face.m_nb[2] = -1;
    face.m_alive = true;
    face.m_visible = false;
    face.m_visit = 0;

    const cVector3d& p0 = a_points[a_v0];
    const cVector3d& p1 = a_points[a_v1];
    const cVector3d& p2 = a_points[a_v2];
    face.m_normal = cNormalize(cCross(p1 - p0, p2 - p0));
    face.m_offset = cDot(face.m_normal, (p0 + p1 + p2) / 3.0);

    return (face);
}


//! Assigns a point to the outside set of the face it is the farthest in front of.
void assignPoint(std::vector<cHullFace>& a_faces,
                 const std::vector<int>& a_candidates,
                 const cVector3d& a_point,
                 int a_index,
                 double a_tolerance)
{
    int best = -1;
    double bestDist = a_tolerance;
    for (size_t i=0; i<a_candidates.size(); i++)
    {
        double d = a_faces[a_candidates[i]].distance(a_point);
        if (d > bestDist)
        {
            bestDist = d;
            best = a_candidates[i];
        }
    }

    // points that are not in front of any face are inside the hull
    if (best >= 0)
    {
        a_faces[best].m_outside.push_back(a_index);
    }
}

//------------------------------------------------------------------------------
} // anonymous namespace
//------------------------------------------------------------------------------


//==============================================================================
/*!
    This method clears the hull.
*/
//==============================================================================
void cConvexHull::clear()
{
    m_vertices.clear();
    m_triangles.clear();
    m_dimension = -1;
    m_tolerance = 0.0;
}


//==============================================================================
/*!
    This method computes the convex hull of the point set \p a_points.

    \param  a_points  Input points.

    \return __true__ if a hull was computed, __false__ if the input is empty.
*/
//==============================================================================
bool cConvexHull::compute(const std::vector<cVector3d>& a_points)
{
    clear();

    const std::vector<cVector3d>& P = a_points;
    const int n = (int)P.size();
    if (n == 0)
    {
        return (false);
    }

    //--------------------------------------------------------------------------
    // TOLERANCE AND EXTREME POINTS
    //--------------------------------------------------------------------------

    // extreme points along each axis; ties are broken lexicographically on the
    // other coordinates, which makes every selected point a hull vertex
    auto lexLess = [&](int a_a, int a_b, int a_axis)
    {
        for (int j=0; j<3; j++)
        {
            int c = (a_axis + j) % 3;
            if (P[a_a](c) != P[a_b](c)) return (P[a_a](c) < P[a_b](c));
        }
        return (false);
    };

    // round-off error of plane distances grows with the magnitude of the
    // coordinates (same estimate as qhull / J. Lloyd's quickhull3d)
    double maxAbs[3] = { 0.0, 0.0, 0.0 };
    int minIdx[3] = { 0, 0, 0 };
    int maxIdx[3] = { 0, 0, 0 };
    for (int i=0; i<n; i++)
    {
        for (int k=0; k<3; k++)
        {
            maxAbs[k] = cMax(maxAbs[k], fabs(P[i](k)));
            if (lexLess(i, minIdx[k], k)) minIdx[k] = i;
            if (lexLess(maxIdx[k], i, k)) maxIdx[k] = i;
        }
    }
    m_tolerance = 3.0 * DBL_EPSILON * (maxAbs[0] + maxAbs[1] + maxAbs[2]);
    const double tol = m_tolerance;


    //--------------------------------------------------------------------------
    // INITIAL SIMPLEX
    //--------------------------------------------------------------------------

    // two most distant extreme points
    int i0 = 0;
    int i1 = 0;
    double best = -1.0;
    for (int k=0; k<3; k++)
    {
        double d = cDistance(P[minIdx[k]], P[maxIdx[k]]);
        if (d > best)
        {
            best = d;
            i0 = minIdx[k];
            i1 = maxIdx[k];
        }
    }
    if (best <= tol)
    {
        // all points coincide
        m_vertices.push_back(P[i0]);
        m_dimension = 0;
        return (true);
    }

    // point farthest from line (i0, i1); ties: farthest along the line
    const cVector3d dir = cNormalize(P[i1] - P[i0]);
    const cVector3d mid = 0.5 * (P[i0] + P[i1]);
    cFarthest far2;
    for (int i=0; i<n; i++)
    {
        far2.test(i, cCross(P[i] - P[i0], dir).length(), fabs(cDot(P[i] - mid, dir)), tol);
    }
    const int i2 = far2.m_index;
    if (far2.m_dist <= tol)
    {
        computeLinear(P, P[i0], dir);
        return (true);
    }

    // point farthest from plane (i0, i1, i2); ties: farthest from i0
    const cVector3d normal = cNormalize(cCross(P[i1] - P[i0], P[i2] - P[i0]));
    cFarthest far3;
    for (int i=0; i<n; i++)
    {
        far3.test(i, fabs(cDot(P[i] - P[i0], normal)), cDistance(P[i], P[i0]), tol);
    }
    const int i3 = far3.m_index;
    if (far3.m_dist <= tol)
    {
        computePlanar(P, P[i0], P[i1], P[i2]);
        return (true);
    }

    // tetrahedron: each face leaves out one vertex, which must lie behind it
    std::vector<cHullFace> faces;
    faces.reserve(8 * 64);
    const int tet[4] = { i0, i1, i2, i3 };
    for (int f=0; f<4; f++)
    {
        int v[3];
        int m = 0;
        for (int j=0; j<4; j++)
        {
            if (j != f) v[m++] = tet[j];
        }
        cHullFace face = newFace(P, v[0], v[1], v[2]);
        if (face.distance(P[tet[f]]) > 0.0)
        {
            face = newFace(P, v[0], v[2], v[1]);
        }
        faces.push_back(face);
    }

    // adjacency of the tetrahedron
    for (int f=0; f<4; f++)
    {
        for (int k=0; k<3; k++)
        {
            int a = faces[f].m_v[k];
            int b = faces[f].m_v[(k+1)%3];
            for (int g=0; g<4; g++)
            {
                for (int j=0; j<3; j++)
                {
                    if ((faces[g].m_v[j] == b) && (faces[g].m_v[(j+1)%3] == a))
                    {
                        faces[f].m_nb[k] = g;
                    }
                }
            }
        }
    }

    // distribute the remaining points among the outside sets
    std::vector<int> candidates;
    for (int f=0; f<4; f++) candidates.push_back(f);
    for (int i=0; i<n; i++)
    {
        if ((i == i0) || (i == i1) || (i == i2) || (i == i3)) continue;
        assignPoint(faces, candidates, P[i], i, tol);
    }


    //--------------------------------------------------------------------------
    // QUICKHULL ITERATIONS
    //--------------------------------------------------------------------------

    // horizon links, indexed by point: -1 when unused
    std::vector<int> startLink(n, -1);
    std::vector<int> endLink(n, -1);

    std::vector<int> stack;
    std::vector<int> visibleFaces;
    std::vector<int> orphans;
    std::vector<cHorizonEdge> horizon;
    unsigned int visitStamp = 0;

    // faces are only appended, and a face receives points only when it is
    // created, so a single forward pass processes every outside set
    for (size_t fi=0; fi<faces.size(); fi++)
    {
        while (faces[fi].m_alive && !faces[fi].m_outside.empty())
        {
            // eye point: the farthest point in front of the face; ties: the
            // farthest from the face centroid
            std::vector<int>& outside = faces[fi].m_outside;
            const cVector3d centroid = (P[faces[fi].m_v[0]] + P[faces[fi].m_v[1]] + P[faces[fi].m_v[2]]) / 3.0;
            cFarthest farEye;
            for (size_t s=0; s<outside.size(); s++)
            {
                const cVector3d& p = P[outside[s]];
                farEye.test((int)s, faces[fi].distance(p), cDistance(p, centroid), tol);
            }
            const size_t eyeSlot = (size_t)farEye.m_index;
            const int eye = outside[eyeSlot];
            const cVector3d& eyePos = P[eye];

            // flood fill the faces visible from the eye and collect the horizon
            visitStamp++;
            visibleFaces.clear();
            horizon.clear();
            faces[fi].m_visit = visitStamp;
            faces[fi].m_visible = true;
            stack.assign(1, (int)fi);
            while (!stack.empty())
            {
                int f = stack.back();
                stack.pop_back();
                visibleFaces.push_back(f);

                for (int k=0; k<3; k++)
                {
                    int g = faces[f].m_nb[k];
                    if (faces[g].m_visit != visitStamp)
                    {
                        faces[g].m_visit = visitStamp;
                        faces[g].m_visible = (faces[g].distance(eyePos) > tol);
                        if (faces[g].m_visible)
                        {
                            stack.push_back(g);
                        }
                    }
                    if (!faces[g].m_visible)
                    {
                        cHorizonEdge edge;
                        edge.m_a = faces[f].m_v[k];
                        edge.m_b = faces[f].m_v[(k+1)%3];
                        edge.m_face = g;
                        horizon.push_back(edge);
                    }
                }
            }

            // the horizon must be one simple closed loop; otherwise the eye
            // point is within round-off of the hull and is dropped
            bool simple = !horizon.empty();
            for (size_t h=0; h<horizon.size(); h++)
            {
                if (startLink[horizon[h].m_a] != -1) simple = false;
                startLink[horizon[h].m_a] = (int)h;
            }
            if (simple)
            {
                size_t count = 0;
                int h = 0;
                do
                {
                    h = startLink[horizon[h].m_b];
                    count++;
                }
                while ((h > 0) && (count <= horizon.size()));
                simple = ((h == 0) && (count == horizon.size()));
            }
            for (size_t h=0; h<horizon.size(); h++)
            {
                startLink[horizon[h].m_a] = -1;
            }
            if (!simple)
            {
                outside[eyeSlot] = outside.back();
                outside.pop_back();
                continue;
            }

            // remove the visible faces and keep their points
            orphans.clear();
            for (size_t i=0; i<visibleFaces.size(); i++)
            {
                cHullFace& face = faces[visibleFaces[i]];
                face.m_alive = false;
                for (size_t s=0; s<face.m_outside.size(); s++)
                {
                    if (face.m_outside[s] != eye) orphans.push_back(face.m_outside[s]);
                }
                std::vector<int>().swap(face.m_outside);
            }

            // cone of new faces from the horizon to the eye point
            const int base = (int)faces.size();
            candidates.clear();
            for (size_t h=0; h<horizon.size(); h++)
            {
                const cHorizonEdge& edge = horizon[h];
                const int id = base + (int)h;

                cHullFace face = newFace(P, edge.m_a, edge.m_b, eye);
                face.m_nb[0] = edge.m_face;
                faces.push_back(face);
                candidates.push_back(id);

                // the non-visible face now borders the new face
                cHullFace& other = faces[edge.m_face];
                for (int k=0; k<3; k++)
                {
                    if ((other.m_v[k] == edge.m_b) && (other.m_v[(k+1)%3] == edge.m_a))
                    {
                        other.m_nb[k] = id;
                    }
                }

                startLink[edge.m_a] = id;
                endLink[edge.m_b] = id;
            }

            // new face (a, b, eye): edge (b, eye) borders the new face starting
            // at b, edge (eye, a) borders the new face ending at a
            for (size_t h=0; h<horizon.size(); h++)
            {
                cHullFace& face = faces[base + h];
                face.m_nb[1] = startLink[face.m_v[1]];
                face.m_nb[2] = endLink[face.m_v[0]];
            }
            for (size_t h=0; h<horizon.size(); h++)
            {
                startLink[horizon[h].m_a] = -1;
                endLink[horizon[h].m_b] = -1;
            }

            // redistribute the points of the removed faces
            for (size_t i=0; i<orphans.size(); i++)
            {
                assignPoint(faces, candidates, P[orphans[i]], orphans[i], tol);
            }
        }
    }


    //--------------------------------------------------------------------------
    // OUTPUT
    //--------------------------------------------------------------------------

    std::vector<int> remap(n, -1);
    for (size_t f=0; f<faces.size(); f++)
    {
        if (!faces[f].m_alive) continue;
        for (int k=0; k<3; k++)
        {
            int index = faces[f].m_v[k];
            if (remap[index] < 0)
            {
                remap[index] = (int)m_vertices.size();
                m_vertices.push_back(P[index]);
            }
            m_triangles.push_back((unsigned int)remap[index]);
        }
    }
    m_dimension = 3;

    return (true);
}


//==============================================================================
/*!
    This method builds the hull of a coplanar point set: a convex polygon
    (Andrew's monotone chain in the plane) returned as a triangle fan.

    \param  a_points  Input points.
    \param  a_p0      First point defining the plane.
    \param  a_p1      Second point defining the plane.
    \param  a_p2      Third point defining the plane.
*/
//==============================================================================
void cConvexHull::computePlanar(const std::vector<cVector3d>& a_points,
                                const cVector3d& a_p0,
                                const cVector3d& a_p1,
                                const cVector3d& a_p2)
{
    const int n = (int)a_points.size();

    // 2D frame of the plane
    cVector3d normal = cNormalize(cCross(a_p1 - a_p0, a_p2 - a_p0));
    cVector3d u = cNormalize(a_p1 - a_p0);
    cVector3d v = cCross(normal, u);

    std::vector<double> x(n), y(n);
    std::vector<int> order(n);
    for (int i=0; i<n; i++)
    {
        cVector3d d = a_points[i] - a_p0;
        x[i] = cDot(d, u);
        y[i] = cDot(d, v);
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b)
    {
        return ((x[a] < x[b]) || ((x[a] == x[b]) && (y[a] < y[b])));
    });

    // __true__ if b is not a strict left turn from a to c: b lies to the
    // right of chord (a, c), or within the tolerance of it
    const double tol = m_tolerance;
    auto notConvex = [&](int a, int b, int c)
    {
        double cross = (x[b] - x[a]) * (y[c] - y[a]) - (y[b] - y[a]) * (x[c] - x[a]);
        return (cross <= tol * sqrt((x[c] - x[a]) * (x[c] - x[a]) + (y[c] - y[a]) * (y[c] - y[a])));
    };

    // monotone chain: lower hull then upper hull, counter clockwise around normal
    std::vector<int> chain(2 * n);
    int k = 0;
    for (int i=0; i<n; i++)
    {
        while ((k >= 2) && notConvex(chain[k-2], chain[k-1], order[i])) k--;
        chain[k++] = order[i];
    }
    for (int i=n-2, lower=k+1; i>=0; i--)
    {
        while ((k >= lower) && notConvex(chain[k-2], chain[k-1], order[i])) k--;
        chain[k++] = order[i];
    }
    chain.resize(cMax(k - 1, 0));

    if (chain.size() < 3)
    {
        computeLinear(a_points, a_p0, u);
        return;
    }

    for (size_t i=0; i<chain.size(); i++)
    {
        m_vertices.push_back(a_points[chain[i]]);
    }
    for (unsigned int i=1; i+1<(unsigned int)chain.size(); i++)
    {
        m_triangles.push_back(0);
        m_triangles.push_back(i);
        m_triangles.push_back(i + 1);
    }
    m_dimension = 2;
}


//==============================================================================
/*!
    This method builds the hull of a collinear point set: the segment joining
    the two extreme points along the line.

    \param  a_points  Input points.
    \param  a_p0      Point on the line.
    \param  a_dir     Unit direction of the line.
*/
//==============================================================================
void cConvexHull::computeLinear(const std::vector<cVector3d>& a_points,
                                const cVector3d& a_p0,
                                const cVector3d& a_dir)
{
    size_t iMin = 0;
    size_t iMax = 0;
    double tMin = cDot(a_points[0] - a_p0, a_dir);
    double tMax = tMin;
    for (size_t i=1; i<a_points.size(); i++)
    {
        double t = cDot(a_points[i] - a_p0, a_dir);
        if (t < tMin) { tMin = t; iMin = i; }
        if (t > tMax) { tMax = t; iMax = i; }
    }

    m_vertices.push_back(a_points[iMin]);
    m_vertices.push_back(a_points[iMax]);
    m_dimension = 1;
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
