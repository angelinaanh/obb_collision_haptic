//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBIntersect.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CCollisionOBBIntersect.h"
//------------------------------------------------------------------------------
#include "math/CMaths.h"
//------------------------------------------------------------------------------
#include <cmath>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//==============================================================================
/*!
    Separating axis test between two oriented boxes. All quantities are
    expressed in the frame of box A: R = A^T B, t = A^T (cB - cA).
*/
//==============================================================================
bool cTestOBBOBB(const cOBBShape& a_boxA, const cOBBShape& a_boxB)
{
    // guards against false separations when two edges are (nearly) parallel
    const double EPS = 1e-6;

    const cVector3d a[3] = { a_boxA.m_rot.getCol0(), a_boxA.m_rot.getCol1(), a_boxA.m_rot.getCol2() };
    const cVector3d b[3] = { a_boxB.m_rot.getCol0(), a_boxB.m_rot.getCol1(), a_boxB.m_rot.getCol2() };
    const double ha[3] = { a_boxA.m_half(0), a_boxA.m_half(1), a_boxA.m_half(2) };
    const double hb[3] = { a_boxB.m_half(0), a_boxB.m_half(1), a_boxB.m_half(2) };

    double R[3][3];
    double AR[3][3];
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            R[i][j] = cDot(a[i], b[j]);
            AR[i][j] = fabs(R[i][j]) + EPS;
        }
    }

    const cVector3d d = a_boxB.m_center - a_boxA.m_center;
    const double t[3] = { cDot(d, a[0]), cDot(d, a[1]), cDot(d, a[2]) };

    double ra, rb;

    // axes of A
    for (int i=0; i<3; i++)
    {
        ra = ha[i];
        rb = hb[0] * AR[i][0] + hb[1] * AR[i][1] + hb[2] * AR[i][2];
        if (fabs(t[i]) > ra + rb) return (false);
    }

    // axes of B
    for (int j=0; j<3; j++)
    {
        ra = ha[0] * AR[0][j] + ha[1] * AR[1][j] + ha[2] * AR[2][j];
        rb = hb[j];
        if (fabs(t[0] * R[0][j] + t[1] * R[1][j] + t[2] * R[2][j]) > ra + rb) return (false);
    }

    // A0 x B0, A0 x B1, A0 x B2
    ra = ha[1] * AR[2][0] + ha[2] * AR[1][0];
    rb = hb[1] * AR[0][2] + hb[2] * AR[0][1];
    if (fabs(t[2] * R[1][0] - t[1] * R[2][0]) > ra + rb) return (false);

    ra = ha[1] * AR[2][1] + ha[2] * AR[1][1];
    rb = hb[0] * AR[0][2] + hb[2] * AR[0][0];
    if (fabs(t[2] * R[1][1] - t[1] * R[2][1]) > ra + rb) return (false);

    ra = ha[1] * AR[2][2] + ha[2] * AR[1][2];
    rb = hb[0] * AR[0][1] + hb[1] * AR[0][0];
    if (fabs(t[2] * R[1][2] - t[1] * R[2][2]) > ra + rb) return (false);

    // A1 x B0, A1 x B1, A1 x B2
    ra = ha[0] * AR[2][0] + ha[2] * AR[0][0];
    rb = hb[1] * AR[1][2] + hb[2] * AR[1][1];
    if (fabs(t[0] * R[2][0] - t[2] * R[0][0]) > ra + rb) return (false);

    ra = ha[0] * AR[2][1] + ha[2] * AR[0][1];
    rb = hb[0] * AR[1][2] + hb[2] * AR[1][0];
    if (fabs(t[0] * R[2][1] - t[2] * R[0][1]) > ra + rb) return (false);

    ra = ha[0] * AR[2][2] + ha[2] * AR[0][2];
    rb = hb[0] * AR[1][1] + hb[1] * AR[1][0];
    if (fabs(t[0] * R[2][2] - t[2] * R[0][2]) > ra + rb) return (false);

    // A2 x B0, A2 x B1, A2 x B2
    ra = ha[0] * AR[1][0] + ha[1] * AR[0][0];
    rb = hb[1] * AR[2][2] + hb[2] * AR[2][1];
    if (fabs(t[1] * R[0][0] - t[0] * R[1][0]) > ra + rb) return (false);

    ra = ha[0] * AR[1][1] + ha[1] * AR[0][1];
    rb = hb[0] * AR[2][2] + hb[2] * AR[2][0];
    if (fabs(t[1] * R[0][1] - t[0] * R[1][1]) > ra + rb) return (false);

    ra = ha[0] * AR[1][2] + ha[1] * AR[0][2];
    rb = hb[0] * AR[2][1] + hb[1] * AR[2][0];
    if (fabs(t[1] * R[0][2] - t[0] * R[1][2]) > ra + rb) return (false);

    // no separating axis found
    return (true);
}


//==============================================================================
/*!
    Slab test: the segment is clipped against the three pairs of parallel
    planes of the box, in the frame of the box.
*/
//==============================================================================
bool cTestSegmentOBB(const cVector3d& a_pointA,
                     const cVector3d& a_pointB,
                     const cOBBShape& a_box,
                     double a_inflate,
                     double& a_tEnter)
{
    const cVector3d p = a_pointA - a_box.m_center;
    const cVector3d d = a_pointB - a_pointA;

    const cVector3d axes[3] = { a_box.m_rot.getCol0(), a_box.m_rot.getCol1(), a_box.m_rot.getCol2() };

    double tMin = 0.0;
    double tMax = 1.0;
    for (int k=0; k<3; k++)
    {
        const cVector3d& axis = axes[k];
        const double pk = cDot(p, axis);
        const double dk = cDot(d, axis);
        const double h = a_box.m_half(k) + a_inflate;

        if (dk == 0.0)
        {
            // segment parallel to the slab
            if (fabs(pk) > h) return (false);
        }
        else
        {
            double t1 = (-h - pk) / dk;
            double t2 = ( h - pk) / dk;
            if (t1 > t2) { double tmp = t1; t1 = t2; t2 = tmp; }
            if (t1 > tMin) tMin = t1;
            if (t2 < tMax) tMax = t2;
            if (tMin > tMax) return (false);
        }
    }

    a_tEnter = tMin;
    return (true);
}


//==============================================================================
/*!
    Moller-Trumbore segment-triangle intersection, accepting both sides.
*/
//==============================================================================
bool cTestSegmentTriangle(const cVector3d& a_pointA,
                          const cVector3d& a_pointB,
                          const cVector3d& a_v0,
                          const cVector3d& a_v1,
                          const cVector3d& a_v2,
                          double& a_t)
{
    const cVector3d dir = a_pointB - a_pointA;
    const cVector3d e1 = a_v1 - a_v0;
    const cVector3d e2 = a_v2 - a_v0;

    const cVector3d pvec = cCross(dir, e2);
    const double det = cDot(e1, pvec);
    if (det == 0.0)
    {
        return (false);
    }
    const double invDet = 1.0 / det;

    const cVector3d tvec = a_pointA - a_v0;
    const double u = cDot(tvec, pvec) * invDet;
    if ((u < 0.0) || (u > 1.0)) return (false);

    const cVector3d qvec = cCross(tvec, e1);
    const double v = cDot(dir, qvec) * invDet;
    if ((v < 0.0) || (u + v > 1.0)) return (false);

    const double t = cDot(e2, qvec) * invDet;
    if ((t < 0.0) || (t > 1.0)) return (false);

    a_t = t;
    return (true);
}


//==============================================================================
/*!
    Separating axis test between two triangles.
*/
//==============================================================================
bool cTestTriangleTriangle(const cVector3d& a_p0,
                           const cVector3d& a_p1,
                           const cVector3d& a_p2,
                           const cVector3d& a_q0,
                           const cVector3d& a_q1,
                           const cVector3d& a_q2)
{
    const cVector3d P[3] = { a_p0, a_p1, a_p2 };
    const cVector3d Q[3] = { a_q0, a_q1, a_q2 };
    const cVector3d ep[3] = { a_p1 - a_p0, a_p2 - a_p1, a_p0 - a_p2 };
    const cVector3d eq[3] = { a_q1 - a_q0, a_q2 - a_q1, a_q0 - a_q2 };
    const cVector3d np = cCross(ep[0], ep[1]);
    const cVector3d nq = cCross(eq[0], eq[1]);

    // __true__ if the projections of the triangles on the axis are disjoint
    auto separated = [&](const cVector3d& a_axis)
    {
        double pMin = cDot(P[0], a_axis), pMax = pMin;
        double qMin = cDot(Q[0], a_axis), qMax = qMin;
        for (int i=1; i<3; i++)
        {
            const double tp = cDot(P[i], a_axis);
            const double tq = cDot(Q[i], a_axis);
            pMin = cMin(pMin, tp); pMax = cMax(pMax, tp);
            qMin = cMin(qMin, tq); qMax = cMax(qMax, tq);
        }
        return ((pMax < qMin) || (qMax < pMin));
    };

    // face normals
    if (separated(np)) return (false);
    if (separated(nq)) return (false);

    // edge-edge axes; a (near) zero cross product carries no direction
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
        {
            const cVector3d axis = cCross(ep[i], eq[j]);
            if (axis.lengthsq() > 1e-24 * ep[i].lengthsq() * eq[j].lengthsq())
            {
                if (separated(axis)) return (false);
            }
        }
    }

    // in-plane edge normals (coplanar triangles, parallel edges)
    for (int i=0; i<3; i++)
    {
        if (separated(cCross(np, ep[i]))) return (false);
        if (separated(cCross(nq, eq[i]))) return (false);
    }

    return (true);
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
