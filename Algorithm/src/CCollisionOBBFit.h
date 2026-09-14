//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBFit.h
    \ingroup    obb

    \brief
    Computes a tight oriented bounding box with convex hull preprocessing and
    continuous (area-weighted) PCA.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBFitH
#define CCollisionOBBFitH
//------------------------------------------------------------------------------
#include "CCollisionOBBBox.h"
#include "CCollisionOBBRefine.h"
#include "CConvexHull.h"
//------------------------------------------------------------------------------
#include <vector>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

class cMesh;

//==============================================================================
/*!
    \file       CCollisionOBBFit.h

    \details
    Fitting pipeline implemented by \ref cComputeOBB(): \n

    1. __Convex hull preprocessing__: only the surface of the convex hull of
       the points is used (\ref cConvexHull). Interior vertices, concavities
       and the extra vertices of finely tessellated flat regions are removed. \n

    2. __Continuous PCA__ (Gottschalk, Lin, Manocha - "OBBTree", RAPID): the
       hull triangles are treated as a surface with uniform area density, so
       the weight of each triangle is its area \f$A^i\f$, not its number of
       vertices. This removes the tessellation bias of vertex-based PCA. \n
       Area-weighted centroid:
       \f[ \mu = \frac{1}{A_H} \sum_i A^i c^i, \qquad A_H = \sum_i A^i \f]
       Covariance integrated over the triangles \f$(p^i, q^i, r^i)\f$:
       \f[ C_{jk} = \frac{1}{A_H} \sum_i \frac{A^i}{12}
           \left( 9 c^i_j c^i_k + p^i_j p^i_k + q^i_j q^i_k + r^i_j r^i_k \right)
           - \mu_j \mu_k \f]

    3. __Axes__: the eigenvectors \f$u_0, u_1, u_2\f$ of \f$C\f$, sorted by
       decreasing eigenvalue and arranged as a right-handed frame. \n

    4. __Iterative refinement__ (\ref cRefineOBBAxes()): small rotations of
       the frame about \f$u_0, u_1, u_2\f$ are kept when they reduce the
       volume of the box around the hull vertices. \n

    5. __Extents__: all input points are projected onto the axes; the box is
       the tightest box around the projections for these axes.
*/
//==============================================================================

//------------------------------------------------------------------------------
/*!
    \brief
    Computes the area-weighted centroid and covariance of a triangle surface.

    \param  a_vertices    Vertices of the surface.
    \param  a_triangles   Triangles (3 vertex indices per triangle).
    \param  a_centroid    Output: area-weighted centroid \f$\mu\f$.
    \param  a_covariance  Output: covariance matrix \f$C\f$.
    \param  a_totalArea   Optional output: total area \f$A_H\f$.

    \return __false__ if the surface has no area (outputs are then set to zero).
*/
//------------------------------------------------------------------------------
bool cComputeSurfaceCovariance(const std::vector<cVector3d>& a_vertices,
                               const std::vector<unsigned int>& a_triangles,
                               cVector3d& a_centroid,
                               cMatrix3d& a_covariance,
                               double* a_totalArea = NULL);

//------------------------------------------------------------------------------
/*!
    \brief
    Computes the mean and covariance of a point set (vertex-based PCA).

    \details
    Used as a fallback when the geometry has no area (points or segments),
    and as a reference to compare with the continuous PCA.
*/
//------------------------------------------------------------------------------
void cComputePointCovariance(const std::vector<cVector3d>& a_points,
                             cVector3d& a_mean,
                             cMatrix3d& a_covariance);

//------------------------------------------------------------------------------
/*!
    \brief
    Computes the principal axes of a covariance matrix.

    \param  a_covariance   Symmetric covariance matrix.
    \param  a_eigenvalues  Optional output: eigenvalues, in decreasing order.

    \return Rotation matrix whose columns are the unit eigenvectors
            \f$u_0, u_1, u_2\f$ sorted by decreasing eigenvalue
            (right-handed: \f$u_2 = u_0 \times u_1\f$).
*/
//------------------------------------------------------------------------------
cMatrix3d cComputePrincipalAxes(const cMatrix3d& a_covariance,
                                cVector3d* a_eigenvalues = NULL);

//------------------------------------------------------------------------------
/*!
    \brief
    Computes the tightest box with given axes that contains a point set.

    \details
    The orientation and the center are stored in float precision by
    cCollisionOBBBox. The half-extents are measured from the stored center
    along the stored axes and rounded up, so that every input point passes
    cCollisionOBBBox::contains().

    \param  a_points  Points to enclose.
    \param  a_axes    Rotation matrix whose columns are the box axes.
*/
//------------------------------------------------------------------------------
cCollisionOBBBox cComputeOBBFromAxes(const std::vector<cVector3d>& a_points,
                                     const cMatrix3d& a_axes);

//------------------------------------------------------------------------------
/*!
    \brief
    Computes a tight OBB of a point set (convex hull + continuous PCA +
    iterative refinement).

    \param  a_points  Points to enclose.
    \param  a_hull    Optional output: convex hull used for the fit.
    \param  a_refine  Parameters of the iterative refinement (enabled by default).
    \param  a_stats   Optional output: statistics of the refinement.
*/
//------------------------------------------------------------------------------
cCollisionOBBBox cComputeOBB(const std::vector<cVector3d>& a_points,
                             cConvexHull* a_hull = NULL,
                             const cOBBRefineSettings& a_refine = cOBBRefineSettings(),
                             cOBBRefineStats* a_stats = NULL);

//------------------------------------------------------------------------------
/*!
    \brief
    Computes a tight OBB of a CHAI3D mesh (convex hull + continuous PCA +
    iterative refinement).

    \details
    The vertices used by the allocated triangles of the mesh are enclosed
    (all vertices if the mesh has no triangle). The box is expressed in the
    local reference frame of the mesh.

    \param  a_mesh    Mesh to enclose.
    \param  a_hull    Optional output: convex hull used for the fit.
    \param  a_refine  Parameters of the iterative refinement (enabled by default).
    \param  a_stats   Optional output: statistics of the refinement.
*/
//------------------------------------------------------------------------------
cCollisionOBBBox cComputeOBB(cMesh* a_mesh,
                             cConvexHull* a_hull = NULL,
                             const cOBBRefineSettings& a_refine = cOBBRefineSettings(),
                             cOBBRefineStats* a_stats = NULL);

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
