//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBRefine.h
    \ingroup    obb

    \brief
    Iterative refinement of OBB axes toward the minimum volume box.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBRefineH
#define CCollisionOBBRefineH
//------------------------------------------------------------------------------
#include "math/CMatrix3d.h"
#include "math/CVector3d.h"
//------------------------------------------------------------------------------
#include <vector>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//==============================================================================
/*!
    \struct     cOBBRefineSettings
    \ingroup    obb

    \brief
    Parameters of the iterative refinement of the OBB axes.

    \details
    Each iteration visits the axes \f$u_0, u_1, u_2\f$ in turn. For every
    angle step \f$\Delta\theta\f$ of \ref m_anglesDeg (coarse to fine), the
    frame is rotated about the current axis by \f$+\Delta\theta\f$ and
    \f$-\Delta\theta\f$; the better rotation is kept if it reduces the volume
    of the box. The refinement stops after \ref m_maxIterations iterations,
    or earlier when a whole iteration brings no improvement.
*/
//==============================================================================
struct cOBBRefineSettings
{
    //! Constructor of cOBBRefineSettings (default: 5 iterations, steps of 5 and 1 degrees).
    cOBBRefineSettings() : m_enabled(true), m_maxIterations(5)
    {
        m_anglesDeg.push_back(5.0);
        m_anglesDeg.push_back(1.0);
    }

    //! If __false__, the PCA axes are used without refinement.
    bool m_enabled;

    //! Maximum number of iterations (one iteration visits the 3 axes).
    int m_maxIterations;

    //! Rotation steps in degrees, tried in this order (coarse to fine).
    std::vector<double> m_anglesDeg;
};


//==============================================================================
/*!
    \struct     cOBBRefineStats
    \ingroup    obb

    \brief
    Statistics reported by the iterative refinement.
*/
//==============================================================================
struct cOBBRefineStats
{
    //! Constructor of cOBBRefineStats.
    cOBBRefineStats() : m_iterations(0), m_acceptedRotations(0), m_evaluations(0),
                        m_initialVolume(0.0), m_finalVolume(0.0) {}

    //! Number of iterations performed.
    int m_iterations;

    //! Number of rotations that were accepted.
    int m_acceptedRotations;

    //! Number of candidate frames evaluated.
    int m_evaluations;

    //! Volume of the box with the initial (PCA) axes.
    double m_initialVolume;

    //! Volume of the box with the refined axes.
    double m_finalVolume;
};


//------------------------------------------------------------------------------
/*!
    \brief
    Computes the size (width, height, depth) of the tightest box with given
    axes around a point set.

    \param  a_points  Points to enclose.
    \param  a_axes    Rotation matrix whose columns are the box axes.

    \return Full lengths of the box along each axis.
*/
//------------------------------------------------------------------------------
cVector3d cComputeBoxSize(const std::vector<cVector3d>& a_points,
                          const cMatrix3d& a_axes);

//------------------------------------------------------------------------------
/*!
    \brief
    Refines OBB axes by small rotations that reduce the volume of the box.

    \details
    Only the extremal points need to be passed (typically the vertices of the
    convex hull): the extents of the box around a point set are those of its
    convex hull. \n

    A candidate frame is accepted when it reduces the volume. When the volumes
    are equal within round-off (flat point sets, whose volume is always zero),
    the candidate is accepted if it reduces the surface area, so that the
    bounding rectangle of flat geometry is refined as well.

    \param  a_points    Extremal points (e.g. convex hull vertices).
    \param  a_axes      Initial axes (columns of a rotation matrix), e.g. from PCA.
    \param  a_settings  Refinement parameters.
    \param  a_stats     Optional output: statistics of the refinement.

    \return Refined axes, as the columns of a right-handed rotation matrix.
*/
//------------------------------------------------------------------------------
cMatrix3d cRefineOBBAxes(const std::vector<cVector3d>& a_points,
                         const cMatrix3d& a_axes,
                         const cOBBRefineSettings& a_settings = cOBBRefineSettings(),
                         cOBBRefineStats* a_stats = NULL);

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
