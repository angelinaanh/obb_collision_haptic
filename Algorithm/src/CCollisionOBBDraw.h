//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBDraw.h
    \ingroup    obb

    \brief
    Wireframe drawing of oriented bounding boxes with OpenGL.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBDrawH
#define CCollisionOBBDrawH
//------------------------------------------------------------------------------
#include "CCollisionOBBBox.h"
#include "graphics/CColor.h"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//! Corner pairs of the 12 edges of a box (corners numbered as in cGetOBBCorners()).
extern const int C_OBB_EDGES[12][2];

//------------------------------------------------------------------------------
/*!
    \brief
    Computes the 8 corners of an oriented box.

    \details
    Corner \p i is \f$c + s_0 h_0 u_0 + s_1 h_1 u_1 + s_2 h_2 u_2\f$ where
    \f$s_k = +1\f$ if bit \p k of \p i is set and \f$-1\f$ otherwise; \f$c\f$,
    \f$h_k\f$ and \f$u_k\f$ are the center, half-extents and axes of the box
    (the columns of cCollisionOBBBox::getRotationMatrix()).
*/
//------------------------------------------------------------------------------
void cGetOBBCorners(const cCollisionOBBBox& a_box, cVector3d a_corners[8]);

//------------------------------------------------------------------------------
/*!
    \brief
    Draws the 12 edges of an oriented box with the current OpenGL color, in
    the current OpenGL frame.
*/
//------------------------------------------------------------------------------
void cDrawWireOBB(const cCollisionOBBBox& a_box);

//------------------------------------------------------------------------------
/*!
    \brief
    Color of a tree level for display: blue at the root, then cyan, green,
    yellow and red at the deepest level.

    \param  a_level  Level normalized to [0, 1] (depth / maximum depth).
*/
//------------------------------------------------------------------------------
cColorf cOBBDepthColor(double a_level);

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
