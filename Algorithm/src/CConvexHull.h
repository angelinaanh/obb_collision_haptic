//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CConvexHull.h
    \ingroup    obb

    \brief
    Implements a 3D convex hull (Quickhull) used to preprocess the geometry
    before fitting an OBB.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CConvexHullH
#define CConvexHullH
//------------------------------------------------------------------------------
#include "math/CVector3d.h"
//------------------------------------------------------------------------------
#include <vector>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//==============================================================================
/*!
    \class      cConvexHull
    \ingroup    obb

    \brief
    This class computes the convex hull of a 3D point set.

    \details
    The hull is computed with the Quickhull algorithm (Barber, Dobkin and
    Huhdanpaa, 1996). Only the points that lie on the hull are kept, so all
    interior vertices of a mesh (and all vertices lying inside a flat region,
    for instance the extra vertices of a finely tessellated plane) are
    discarded. \n

    Depending on the input, the hull has one of the following dimensions: \n

    - __3__: a closed triangulated polyhedron. Triangles are oriented counter
      clockwise when seen from outside (outward normals). \n
    - __2__: all points are coplanar. The hull is a convex polygon returned
      as a triangle fan (one side only). \n
    - __1__: all points are collinear. The hull is a segment (2 vertices, no
      triangles). \n
    - __0__: all points coincide (1 vertex, no triangles). \n
    - __-1__: the input is empty. \n

    Two points are considered equal, and a point is considered to lie on a
    plane, when their distance is below a tolerance that scales with the size
    of the coordinates (see \ref getTolerance()).
*/
//==============================================================================
class cConvexHull
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of cConvexHull.
    cConvexHull() { clear(); }

    //! Destructor of cConvexHull.
    virtual ~cConvexHull() {}


    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:

    //! This method computes the convex hull of a point set. Returns __false__ if the input is empty.
    bool compute(const std::vector<cVector3d>& a_points);

    //! This method clears the hull.
    void clear();

    //! This method returns the dimension of the hull (-1, 0, 1, 2 or 3).
    int getDimension() const { return (m_dimension); }

    //! This method returns the vertices of the hull.
    const std::vector<cVector3d>& getVertices() const { return (m_vertices); }

    //! This method returns the triangles of the hull (3 vertex indices per triangle).
    const std::vector<unsigned int>& getTriangles() const { return (m_triangles); }

    //! This method returns the number of vertices of the hull.
    unsigned int getNumVertices() const { return ((unsigned int)m_vertices.size()); }

    //! This method returns the number of triangles of the hull.
    unsigned int getNumTriangles() const { return ((unsigned int)m_triangles.size() / 3); }

    //! This method returns the distance tolerance used during the last computation.
    double getTolerance() const { return (m_tolerance); }


    //--------------------------------------------------------------------------
    // PROTECTED METHODS:
    //--------------------------------------------------------------------------

protected:

    //! This method builds a planar hull (all points lie in the plane through a_p0, a_p1, a_p2).
    void computePlanar(const std::vector<cVector3d>& a_points,
                       const cVector3d& a_p0,
                       const cVector3d& a_p1,
                       const cVector3d& a_p2);

    //! This method builds a segment hull (all points lie on the line through a_p0 along a_dir).
    void computeLinear(const std::vector<cVector3d>& a_points,
                       const cVector3d& a_p0,
                       const cVector3d& a_dir);


    //--------------------------------------------------------------------------
    // PROTECTED MEMBERS:
    //--------------------------------------------------------------------------

protected:

    //! Vertices of the hull.
    std::vector<cVector3d> m_vertices;

    //! Triangles of the hull (3 indices into m_vertices per triangle).
    std::vector<unsigned int> m_triangles;

    //! Dimension of the hull.
    int m_dimension;

    //! Distance tolerance used during the last computation.
    double m_tolerance;
};

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
