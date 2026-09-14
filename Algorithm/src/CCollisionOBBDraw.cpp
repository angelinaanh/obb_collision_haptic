//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBDraw.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CCollisionOBBDraw.h"
//------------------------------------------------------------------------------
#include "math/CMaths.h"
#include "system/CGlobals.h"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

// edges join corners that differ in exactly one bit (one axis)
const int C_OBB_EDGES[12][2] =
{
    {0, 1}, {2, 3}, {4, 5}, {6, 7},     // along u0
    {0, 2}, {1, 3}, {4, 6}, {5, 7},     // along u1
    {0, 4}, {1, 5}, {2, 6}, {3, 7}      // along u2
};


//==============================================================================
/*!
    This function computes the 8 corners of an oriented box.
*/
//==============================================================================
void cGetOBBCorners(const cCollisionOBBBox& a_box, cVector3d a_corners[8])
{
    const cVector3d c = a_box.getCenter();
    const cMatrix3d rot = a_box.getRotationMatrix();
    const cVector3d e0 = (double)a_box.m_halfExtents[0] * rot.getCol0();
    const cVector3d e1 = (double)a_box.m_halfExtents[1] * rot.getCol1();
    const cVector3d e2 = (double)a_box.m_halfExtents[2] * rot.getCol2();

    for (int i=0; i<8; i++)
    {
        a_corners[i] = c + ((i & 1) ? e0 : -e0)
                         + ((i & 2) ? e1 : -e1)
                         + ((i & 4) ? e2 : -e2);
    }
}


//==============================================================================
/*!
    This function draws the edges of an oriented box.
*/
//==============================================================================
void cDrawWireOBB(const cCollisionOBBBox& a_box)
{
#ifdef C_USE_OPENGL

    cVector3d corners[8];
    cGetOBBCorners(a_box, corners);

    glBegin(GL_LINES);
    for (int e=0; e<12; e++)
    {
        const cVector3d& a = corners[C_OBB_EDGES[e][0]];
        const cVector3d& b = corners[C_OBB_EDGES[e][1]];
        glVertex3d(a(0), a(1), a(2));
        glVertex3d(b(0), b(1), b(2));
    }
    glEnd();

#endif
}


//==============================================================================
/*!
    This function returns the display color of a normalized tree level.
*/
//==============================================================================
cColorf cOBBDepthColor(double a_level)
{
    // piecewise linear ramp: blue -> cyan -> green -> yellow -> red
    static const float ramp[5][3] = { {0.2f, 0.4f, 1.0f}, {0.0f, 0.9f, 1.0f}, {0.1f, 1.0f, 0.2f},
                                      {1.0f, 0.9f, 0.0f}, {1.0f, 0.2f, 0.1f} };

    const double t = 4.0 * cClamp(a_level, 0.0, 1.0);
    const int i = cMin((int)t, 3);
    const float f = (float)(t - i);

    return (cColorf(ramp[i][0] + f * (ramp[i+1][0] - ramp[i][0]),
                    ramp[i][1] + f * (ramp[i+1][1] - ramp[i][1]),
                    ramp[i][2] + f * (ramp[i+1][2] - ramp[i][2]),
                    1.0f));
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
