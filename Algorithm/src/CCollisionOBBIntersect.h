//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBIntersect.h
    \ingroup    obb

    \brief
    Intersection tests used by the OBB tree: box-box, segment-box,
    segment-triangle and triangle-triangle.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBIntersectH
#define CCollisionOBBIntersectH
//------------------------------------------------------------------------------
#include "CCollisionOBBBox.h"
#include "math/CMatrix3d.h"
#include "math/CVector3d.h"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//==============================================================================
/*!
    \struct     cOBBShape
    \ingroup    obb

    \brief
    Oriented box unpacked in double precision for intersection tests.

    \details
    cCollisionOBBBox stores its orientation as a float quaternion to stay
    compact; the tests need the axes as a rotation matrix. A cOBBShape is the
    box expressed in a given frame: center, axes (columns of m_rot) and
    half-extents.
*/
//==============================================================================
struct cOBBShape
{
    //! Constructor of cOBBShape (empty box at the origin).
    cOBBShape() : m_center(0, 0, 0), m_rot(cIdentity3d()), m_half(0, 0, 0) {}

    //! Constructor of cOBBShape from a stored box.
    explicit cOBBShape(const cCollisionOBBBox& a_box) :
        m_center(a_box.getCenter()),
        m_rot(a_box.getRotationMatrix()),
        m_half(a_box.getHalfExtents()) {}

    //! This method returns the box moved by the rigid transform x' = a_rot * x + a_pos.
    cOBBShape transformed(const cMatrix3d& a_rot, const cVector3d& a_pos) const
    {
        cOBBShape shape;
        shape.m_center = a_rot * m_center + a_pos;
        shape.m_rot = a_rot * m_rot;
        shape.m_half = m_half;
        return (shape);
    }

    //! Center of the box.
    cVector3d m_center;

    //! Axes of the box (columns).
    cMatrix3d m_rot;

    //! Half-extents of the box along its axes.
    cVector3d m_half;
};


//------------------------------------------------------------------------------
/*!
    \brief
    Tests whether two oriented boxes overlap (separating axis theorem).

    \details
    The 15 candidate axes are the 3 axes of each box and the 9 cross products
    of their axes (Gottschalk, Lin, Manocha - "OBBTree", 1996; implementation
    after C. Ericson, "Real-Time Collision Detection", 4.4.1). A small epsilon
    is added to the rotation terms so that nearly parallel axes can never
    produce a false separation: the test is conservative.

    \return __true__ if the boxes overlap (or touch), __false__ if separated.
*/
//------------------------------------------------------------------------------
bool cTestOBBOBB(const cOBBShape& a_boxA, const cOBBShape& a_boxB);

//------------------------------------------------------------------------------
/*!
    \brief
    Tests whether segment [A, B] intersects an oriented box (slab test).

    \param  a_pointA   Start point of the segment.
    \param  a_pointB   End point of the segment.
    \param  a_box      Box.
    \param  a_inflate  Distance added to every half-extent (e.g. a tool radius).
    \param  a_tEnter   Output: segment parameter in [0, 1] where the segment
                       enters the box (0 if A is inside).

    \return __true__ if the segment intersects the (inflated) box.
*/
//------------------------------------------------------------------------------
bool cTestSegmentOBB(const cVector3d& a_pointA,
                     const cVector3d& a_pointB,
                     const cOBBShape& a_box,
                     double a_inflate,
                     double& a_tEnter);

//------------------------------------------------------------------------------
/*!
    \brief
    Tests whether segment [A, B] intersects a triangle (Moller-Trumbore,
    both sides of the triangle).

    \param  a_t  Output: segment parameter in [0, 1] of the intersection point.

    \return __true__ if the segment crosses the triangle. A segment lying in
            the plane of the triangle is not reported.
*/
//------------------------------------------------------------------------------
bool cTestSegmentTriangle(const cVector3d& a_pointA,
                          const cVector3d& a_pointB,
                          const cVector3d& a_v0,
                          const cVector3d& a_v1,
                          const cVector3d& a_v2,
                          double& a_t);

//------------------------------------------------------------------------------
/*!
    \brief
    Tests whether two triangles intersect (separating axis theorem).

    \details
    Candidate axes: the two triangle normals, the 9 cross products of their
    edges, and the 6 in-plane edge normals (needed when the triangles are
    coplanar or have parallel edges). Touching triangles are reported as
    intersecting.

    \return __true__ if the triangles intersect or touch.
*/
//------------------------------------------------------------------------------
bool cTestTriangleTriangle(const cVector3d& a_p0,
                           const cVector3d& a_p1,
                           const cVector3d& a_p2,
                           const cVector3d& a_q0,
                           const cVector3d& a_q1,
                           const cVector3d& a_q2);

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
