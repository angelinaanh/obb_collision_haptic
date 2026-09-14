//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBBox.h
    \ingroup    obb

    \brief
    Implements the node volume of an oriented bounding box collision tree (OBB).
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBBoxH
#define CCollisionOBBBoxH
//------------------------------------------------------------------------------
#include "math/CMaths.h"
#include "math/CQuaternion.h"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//==============================================================================
/*!
    \class      cCollisionOBBBox
    \ingroup    obb

    \brief
    This structure implements an oriented bounding box (OBB).

    \details
    An OBB is a box of arbitrary orientation. It is described by: \n

    - __Center__: 3 floats (x, y, z), expressed in the reference frame of the
      object that owns the box. \n
    - __Half-extents__: 3 floats, half the size of the box along each of its
      local axes. \n
    - __Rotation__: a unit quaternion stored as 4 floats in the order
      __(x, y, z, w)__. The columns of the rotation matrix of this quaternion
      are the local axes of the box. \n

    The whole box therefore occupies exactly 10 floats (40 bytes). The class
    is intentionally non-virtual so that no vtable pointer is added and large
    arrays of nodes stay compact in memory. \n

    __Note__: cQuaternion of CHAI3D stores its components in the order
    __(w, x, y, z)__ and in double precision. The conversion methods of this
    class take care of both the reordering and the precision change.
*/
//==============================================================================
struct cCollisionOBBBox
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Default constructor of cCollisionOBBBox (zero size, identity rotation).
    cCollisionOBBBox() { setEmpty(); }

    //! Constructor of cCollisionOBBBox (orientation given by a quaternion).
    cCollisionOBBBox(const cVector3d& a_center,
                     const cVector3d& a_halfExtents,
                     const cQuaternion& a_rotation)
    {
        setCenter(a_center);
        setHalfExtents(a_halfExtents);
        setRotation(a_rotation);
    }

    //! Constructor of cCollisionOBBBox (orientation given by a rotation matrix).
    cCollisionOBBBox(const cVector3d& a_center,
                     const cVector3d& a_halfExtents,
                     const cMatrix3d& a_rotation)
    {
        setCenter(a_center);
        setHalfExtents(a_halfExtents);
        setRotation(a_rotation);
    }


    //--------------------------------------------------------------------------
    // METHODS - CENTER:
    //--------------------------------------------------------------------------

public:

    //! This method sets the center of the box.
    inline void setCenter(const float a_x, const float a_y, const float a_z)
    {
        m_center[0] = a_x;
        m_center[1] = a_y;
        m_center[2] = a_z;
    }

    //! This method sets the center of the box.
    inline void setCenter(const cVector3d& a_center)
    {
        setCenter((float)a_center(0), (float)a_center(1), (float)a_center(2));
    }

    //! This method returns the center of the box.
    inline cVector3d getCenter() const
    {
        return (cVector3d(m_center[0], m_center[1], m_center[2]));
    }


    //--------------------------------------------------------------------------
    // METHODS - HALF-EXTENTS:
    //--------------------------------------------------------------------------

public:

    //! This method sets the half-extents (half the size along each local axis) of the box.
    inline void setHalfExtents(const float a_x, const float a_y, const float a_z)
    {
        m_halfExtents[0] = a_x;
        m_halfExtents[1] = a_y;
        m_halfExtents[2] = a_z;
    }

    //! This method sets the half-extents (half the size along each local axis) of the box.
    inline void setHalfExtents(const cVector3d& a_halfExtents)
    {
        setHalfExtents((float)a_halfExtents(0), (float)a_halfExtents(1), (float)a_halfExtents(2));
    }

    //! This method returns the half-extents (half the size along each local axis) of the box.
    inline cVector3d getHalfExtents() const
    {
        return (cVector3d(m_halfExtents[0], m_halfExtents[1], m_halfExtents[2]));
    }


    //--------------------------------------------------------------------------
    // METHODS - ROTATION:
    //--------------------------------------------------------------------------

public:

    //--------------------------------------------------------------------------
    /*!
        \brief
        This method sets the orientation of the box from quaternion components.

        \details
        This method sets the orientation of the box from the quaternion
        components passed in the order __(x, y, z, w)__. The quaternion is
        normalized before being stored.

        \param  a_x  Component __x__ of quaternion.
        \param  a_y  Component __y__ of quaternion.
        \param  a_z  Component __z__ of quaternion.
        \param  a_w  Component __w__ of quaternion.
    */
    //--------------------------------------------------------------------------
    inline void setRotation(const float a_x, const float a_y, const float a_z, const float a_w)
    {
        m_rotation[0] = a_x;
        m_rotation[1] = a_y;
        m_rotation[2] = a_z;
        m_rotation[3] = a_w;
        normalizeRotation();
    }

    //! This method sets the orientation of the box from a CHAI3D quaternion __(w, x, y, z)__.
    inline void setRotation(const cQuaternion& a_rotation)
    {
        setRotation((float)a_rotation.x, (float)a_rotation.y, (float)a_rotation.z, (float)a_rotation.w);
    }

    //! This method sets the orientation of the box from a rotation matrix whose columns are the local axes.
    inline void setRotation(const cMatrix3d& a_rotation)
    {
        cQuaternion q;
        q.fromRotMat(a_rotation);
        setRotation(q);
    }

    //--------------------------------------------------------------------------
    /*!
        \brief
        This method returns the orientation of the box as a CHAI3D quaternion
        __(w, x, y, z)__.

        \details
        The float components are unit only to float precision (~1e-7); the
        quaternion is normalized again in double precision so that the rotation
        matrices derived from it (getRotationMatrix(), getAxis()) are
        orthonormal to double precision.
    */
    //--------------------------------------------------------------------------
    inline cQuaternion getRotation() const
    {
        cQuaternion q(m_rotation[3], m_rotation[0], m_rotation[1], m_rotation[2]);
        const double mag = q.mag();
        if (mag < C_SMALL)
        {
            return (cQuaternion(1.0, 0.0, 0.0, 0.0));
        }
        q.w /= mag;
        q.x /= mag;
        q.y /= mag;
        q.z /= mag;
        return (q);
    }

    //! This method returns the orientation of the box as a rotation matrix whose columns are the local axes.
    inline cMatrix3d getRotationMatrix() const
    {
        cMatrix3d rot;
        getRotation().toRotMat(rot);
        return (rot);
    }


    //--------------------------------------------------------------------------
    /*!
        \brief
        This method returns one of the local axes of the box.

        \details
        This method returns the local axis \p a_index (0 = X, 1 = Y, 2 = Z) of
        the box, expressed in the reference frame of the box's center. The
        returned vector has unit length.

        \param  a_index  Index of the axis (0, 1 or 2).

        \return Local axis of the box.
    */
    //--------------------------------------------------------------------------
    inline cVector3d getAxis(const int a_index) const
    {
        cMatrix3d rot = getRotationMatrix();
        switch (a_index)
        {
            case 0:  return (rot.getCol0());
            case 1:  return (rot.getCol1());
            default: return (rot.getCol2());
        }
    }


    //--------------------------------------------------------------------------
    /*!
        \brief
        This method normalizes the rotation quaternion of the box.

        \details
        This method normalizes the rotation quaternion so that it represents a
        pure rotation. If the quaternion is degenerated (length close to zero),
        it is reset to the identity rotation.
    */
    //--------------------------------------------------------------------------
    inline void normalizeRotation()
    {
        double mag = sqrt((double)m_rotation[0] * m_rotation[0] +
                          (double)m_rotation[1] * m_rotation[1] +
                          (double)m_rotation[2] * m_rotation[2] +
                          (double)m_rotation[3] * m_rotation[3]);

        if (mag < C_SMALL)
        {
            m_rotation[0] = 0.0f;
            m_rotation[1] = 0.0f;
            m_rotation[2] = 0.0f;
            m_rotation[3] = 1.0f;
        }
        else
        {
            for (int i=0; i<4; i++)
            {
                m_rotation[i] = (float)(m_rotation[i] / mag);
            }
        }
    }


    //--------------------------------------------------------------------------
    // METHODS - GENERAL:
    //--------------------------------------------------------------------------

public:

    //! This method resets the box to a zero size box at the origin with identity rotation.
    inline void setEmpty()
    {
        setCenter(0.0f, 0.0f, 0.0f);
        setHalfExtents(0.0f, 0.0f, 0.0f);
        setRotation(0.0f, 0.0f, 0.0f, 1.0f);
    }


    //--------------------------------------------------------------------------
    /*!
        \brief
        This method tests whether this box contains a point passed as argument.

        \details
        This method projects the vector going from the center of the box to
        the point \p a_point onto each local axis and compares the result with
        the corresponding half-extent. Points located on the surface of the
        box are considered inside.

        \param  a_point  Point to be tested.

        \return __true__ if this box contains this point, __false__ otherwise.
    */
    //--------------------------------------------------------------------------
    inline bool contains(const cVector3d& a_point) const
    {
        cVector3d d = a_point - getCenter();
        cMatrix3d rot = getRotationMatrix();

        if (cAbs(cDot(d, rot.getCol0())) > m_halfExtents[0]) return (false);
        if (cAbs(cDot(d, rot.getCol1())) > m_halfExtents[1]) return (false);
        if (cAbs(cDot(d, rot.getCol2())) > m_halfExtents[2]) return (false);

        return (true);
    }


    //--------------------------------------------------------------------------
    // PUBLIC MEMBERS:
    //--------------------------------------------------------------------------

public:

    //! Center of the box (x, y, z).
    float m_center[3];

    //! Half-extents of the box along its local axes X, Y and Z.
    float m_halfExtents[3];

    //! Rotation of the box as a unit quaternion stored in the order (x, y, z, w).
    float m_rotation[4];
};


//------------------------------------------------------------------------------
// Guarantee the compact layout: 3 + 3 + 4 floats, no padding, no vtable.
static_assert(sizeof(cCollisionOBBBox) == 10 * sizeof(float),
              "cCollisionOBBBox must occupy exactly 10 floats.");
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
