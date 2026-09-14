//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBFit.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CCollisionOBBFit.h"
//------------------------------------------------------------------------------
#include "world/CMesh.h"
//------------------------------------------------------------------------------
#include <Eigen/Eigenvalues>
#include <cfloat>
#include <cmath>
#include <limits>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace {
//------------------------------------------------------------------------------

inline Eigen::Vector3d toEigen(const cVector3d& a_vector)
{
    return (Eigen::Vector3d(a_vector(0), a_vector(1), a_vector(2)));
}

//! Mean of a point set, used as reference point to reduce round-off.
Eigen::Vector3d computeMean(const std::vector<cVector3d>& a_points)
{
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (size_t i=0; i<a_points.size(); i++)
    {
        mean += toEigen(a_points[i]);
    }
    return (mean / (double)a_points.size());
}

//! Smallest float that is greater than or equal to a_value.
inline float roundUp(double a_value)
{
    float result = (float)a_value;
    if ((double)result < a_value)
    {
        result = std::nextafter(result, std::numeric_limits<float>::infinity());
    }
    return (result);
}

//------------------------------------------------------------------------------
} // anonymous namespace
//------------------------------------------------------------------------------


//==============================================================================
/*!
    This function computes the area-weighted centroid and the covariance matrix
    of a triangle surface with uniform area density (Gottschalk, RAPID).

    The second moment of a triangle (p, q, r) with centroid c and area A is
    A/12 * (9 c c^T + p p^T + q q^T + r r^T). The coordinates are taken
    relative to the vertex mean, which does not change C but avoids the
    cancellation in E[x x^T] - mu mu^T for surfaces far from the origin.
*/
//==============================================================================
bool cComputeSurfaceCovariance(const std::vector<cVector3d>& a_vertices,
                               const std::vector<unsigned int>& a_triangles,
                               cVector3d& a_centroid,
                               cMatrix3d& a_covariance,
                               double* a_totalArea)
{
    a_centroid.zero();
    a_covariance = cMatrix3d(Eigen::Matrix3d::Zero());
    if (a_totalArea != NULL) *a_totalArea = 0.0;

    if (a_vertices.empty() || (a_triangles.size() < 3))
    {
        return (false);
    }

    const Eigen::Vector3d ref = computeMean(a_vertices);

    double area = 0.0;
    Eigen::Vector3d moment = Eigen::Vector3d::Zero();
    Eigen::Matrix3d second = Eigen::Matrix3d::Zero();

    const size_t numTriangles = a_triangles.size() / 3;
    for (size_t i=0; i<numTriangles; i++)
    {
        const Eigen::Vector3d p = toEigen(a_vertices[a_triangles[3*i+0]]) - ref;
        const Eigen::Vector3d q = toEigen(a_vertices[a_triangles[3*i+1]]) - ref;
        const Eigen::Vector3d r = toEigen(a_vertices[a_triangles[3*i+2]]) - ref;

        const double A = 0.5 * (q - p).cross(r - p).norm();
        const Eigen::Vector3d c = (p + q + r) / 3.0;

        area   += A;
        moment += A * c;
        second += (A / 12.0) * (9.0 * c * c.transpose() +
                                p * p.transpose() +
                                q * q.transpose() +
                                r * r.transpose());
    }

    if (!(area > 0.0))
    {
        return (false);
    }

    const Eigen::Vector3d mu = moment / area;
    const Eigen::Matrix3d C = second / area - mu * mu.transpose();

    a_centroid = cVector3d(mu + ref);
    a_covariance = cMatrix3d(C);
    if (a_totalArea != NULL) *a_totalArea = area;

    return (true);
}


//==============================================================================
/*!
    This function computes the mean and the covariance matrix of a point set.
*/
//==============================================================================
void cComputePointCovariance(const std::vector<cVector3d>& a_points,
                             cVector3d& a_mean,
                             cMatrix3d& a_covariance)
{
    a_mean.zero();
    a_covariance = cMatrix3d(Eigen::Matrix3d::Zero());
    if (a_points.empty())
    {
        return;
    }

    const Eigen::Vector3d mean = computeMean(a_points);
    Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
    for (size_t i=0; i<a_points.size(); i++)
    {
        const Eigen::Vector3d d = toEigen(a_points[i]) - mean;
        C += d * d.transpose();
    }

    a_mean = cVector3d(mean);
    a_covariance = cMatrix3d(C / (double)a_points.size());
}


//==============================================================================
/*!
    This function computes the eigenvectors of a symmetric covariance matrix
    and returns them as the columns of a right-handed rotation matrix, sorted
    by decreasing eigenvalue.
*/
//==============================================================================
cMatrix3d cComputePrincipalAxes(const cMatrix3d& a_covariance,
                                cVector3d* a_eigenvalues)
{
    Eigen::Matrix3d C;
    a_covariance.copyto(C);
    C = 0.5 * (C + C.transpose());

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(C);
    if (solver.info() != Eigen::Success)
    {
        if (a_eigenvalues != NULL) a_eigenvalues->zero();
        return (cIdentity3d());
    }

    // Eigen sorts eigenvalues in increasing order
    const Eigen::Vector3d& values = solver.eigenvalues();
    const Eigen::Matrix3d& vectors = solver.eigenvectors();

    Eigen::Vector3d u0 = vectors.col(2).normalized();
    Eigen::Vector3d u1 = (vectors.col(1) - u0.dot(vectors.col(1)) * u0).normalized();
    Eigen::Vector3d u2 = u0.cross(u1);

    if (a_eigenvalues != NULL)
    {
        a_eigenvalues->set(values(2), values(1), values(0));
    }

    return (cMatrix3d(u0, u1, u2));
}


//==============================================================================
/*!
    This function computes the tightest box with the given axes around a point
    set, taking the float storage of cCollisionOBBBox into account.
*/
//==============================================================================
cCollisionOBBBox cComputeOBBFromAxes(const std::vector<cVector3d>& a_points,
                                     const cMatrix3d& a_axes)
{
    cCollisionOBBBox box;
    if (a_points.empty())
    {
        return (box);
    }

    // use the axes exactly as they are stored (float quaternion)
    box.setRotation(a_axes);
    const cMatrix3d rot = box.getRotationMatrix();
    const cVector3d axis[3] = { rot.getCol0(), rot.getCol1(), rot.getCol2() };

    // range of the projections on each axis, relative to a point of the set:
    // the stored axes are unit only to float precision, and projecting large
    // absolute coordinates would scale that error by the distance to origin
    const cVector3d ref = a_points[0];
    double lo[3] = {  DBL_MAX,  DBL_MAX,  DBL_MAX };
    double hi[3] = { -DBL_MAX, -DBL_MAX, -DBL_MAX };
    for (size_t i=0; i<a_points.size(); i++)
    {
        cVector3d d = a_points[i] - ref;
        for (int k=0; k<3; k++)
        {
            double t = cDot(d, axis[k]);
            lo[k] = cMin(lo[k], t);
            hi[k] = cMax(hi[k], t);
        }
    }

    box.setCenter(ref + 0.5 * (lo[0] + hi[0]) * axis[0] +
                        0.5 * (lo[1] + hi[1]) * axis[1] +
                        0.5 * (lo[2] + hi[2]) * axis[2]);

    // half-extents measured from the stored center, as in contains(), then
    // rounded up to the next float so that no point ends up outside
    const cVector3d center = box.getCenter();
    double h[3] = { 0.0, 0.0, 0.0 };
    for (size_t i=0; i<a_points.size(); i++)
    {
        cVector3d d = a_points[i] - center;
        for (int k=0; k<3; k++)
        {
            h[k] = cMax(h[k], fabs(cDot(d, axis[k])));
        }
    }
    box.setHalfExtents(roundUp(h[0]), roundUp(h[1]), roundUp(h[2]));

    return (box);
}


//==============================================================================
/*!
    This function computes a tight OBB of a point set: convex hull, continuous
    PCA of the hull surface, iterative refinement of the axes, then fit of the
    extents.
*/
//==============================================================================
cCollisionOBBBox cComputeOBB(const std::vector<cVector3d>& a_points,
                             cConvexHull* a_hull,
                             const cOBBRefineSettings& a_refine,
                             cOBBRefineStats* a_stats)
{
    cConvexHull localHull;
    cConvexHull& hull = (a_hull != NULL) ? *a_hull : localHull;

    if (a_stats != NULL) *a_stats = cOBBRefineStats();
    if (!hull.compute(a_points))
    {
        return (cCollisionOBBBox());
    }

    cVector3d centroid;
    cMatrix3d covariance;
    if (!cComputeSurfaceCovariance(hull.getVertices(), hull.getTriangles(), centroid, covariance))
    {
        // hull without area (single point or segment)
        cComputePointCovariance(hull.getVertices(), centroid, covariance);
    }
    cMatrix3d axes = cComputePrincipalAxes(covariance);

    // the extents of a box around the points are those of the hull, so the
    // refinement only needs to project the hull vertices
    axes = cRefineOBBAxes(hull.getVertices(), axes, a_refine, a_stats);

    // extents are fitted on all input points, so the box encloses them even
    // if the hull dropped a point that was within round-off of its surface
    return (cComputeOBBFromAxes(a_points, axes));
}


//==============================================================================
/*!
    This function computes a tight OBB of a mesh, in the local reference frame
    of the mesh.
*/
//==============================================================================
cCollisionOBBBox cComputeOBB(cMesh* a_mesh,
                             cConvexHull* a_hull,
                             const cOBBRefineSettings& a_refine,
                             cOBBRefineStats* a_stats)
{
    if (a_mesh == NULL)
    {
        if (a_stats != NULL) *a_stats = cOBBRefineStats();
        return (cCollisionOBBBox());
    }

    cVertexArrayPtr vertices = a_mesh->m_vertices;
    cTriangleArrayPtr triangles = a_mesh->m_triangles;
    const unsigned int numVertices = vertices->getNumElements();
    const unsigned int numTriangles = triangles->getNumElements();

    // vertices referenced by allocated triangles
    std::vector<cVector3d> points;
    std::vector<bool> used(numVertices, false);
    for (unsigned int i=0; i<numTriangles; i++)
    {
        if (!triangles->getAllocated(i)) continue;

        const unsigned int index[3] = { triangles->getVertexIndex0(i),
                                        triangles->getVertexIndex1(i),
                                        triangles->getVertexIndex2(i) };
        for (int k=0; k<3; k++)
        {
            if (!used[index[k]])
            {
                used[index[k]] = true;
                points.push_back(vertices->getLocalPos(index[k]));
            }
        }
    }

    // mesh without triangles: use all vertices
    if (points.empty())
    {
        for (unsigned int i=0; i<numVertices; i++)
        {
            points.push_back(vertices->getLocalPos(i));
        }
    }

    return (cComputeOBB(points, a_hull, a_refine, a_stats));
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
