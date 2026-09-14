//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBB.h
    \ingroup    obb

    \brief
    CHAI3D collision detector based on the OBB tree, and intersection test
    between two meshes (e.g. a hand and an object).
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBH
#define CCollisionOBBH
//------------------------------------------------------------------------------
#include "CCollisionOBBTree.h"
#include "collisions/CGenericCollision.h"
//------------------------------------------------------------------------------
#include <vector>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

class cMesh;

//==============================================================================
/*!
    \class      cCollisionOBB
    \ingroup    obb

    \brief
    Collision detector for a cMesh based on an OBB tree built with SAH.

    \details
    Drop-in replacement of CHAI3D's AABB detector:

    \code
    cCreateOBBCollisionDetector(mesh);      // builds the tree, attaches it to the mesh
    \endcode

    Segment queries from CHAI3D (haptic tools, cGenericObject::
    computeCollisionDetection(), ...) are culled with the tree; the candidate
    triangles are then tested by the mesh's own triangle array, so the
    collision events (nearest collision, radius of the tool, front/back side
    settings of the material) are exactly those of the other CHAI3D detectors. \n

    The tree also supports tree-tree queries: see cComputeMeshIntersection(). \n

    __Display__. The boxes are drawn as wireframes in the local frame of the
    mesh, so they follow its position and orientation:

    \code
    mesh->setShowCollisionDetector(true);   // show / hide (CHAI3D standard)
    detector->setDisplayDepth(3);           // level 3 of the tree (see below)
    detector->setDisplayDepth(-3);          // all levels from the root to 3
    detector->setColorByDepth(true);        // one color per level
    \endcode

    With a positive display depth d, the nodes of depth d are drawn together
    with the leaves shallower than d: the SAH tree is not uniformly deep, and
    this "cut" through the tree always covers the whole mesh.
*/
//==============================================================================
class cCollisionOBB : public cGenericCollision
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of cCollisionOBB.
    cCollisionOBB();

    //! Destructor of cCollisionOBB.
    virtual ~cCollisionOBB() {}


    //--------------------------------------------------------------------------
    // PUBLIC METHODS:
    //--------------------------------------------------------------------------

public:

    //! This method builds the OBB tree of a mesh.
    bool initialize(cMesh* a_mesh, const cOBBTreeSettings& a_settings = cOBBTreeSettings());

    //! This method rebuilds the tree; it must be called when the mesh is modified.
    virtual void update();

    //! This method computes the collisions between a segment (local coordinates of the mesh) and the mesh.
    virtual bool computeCollision(cGenericObject* a_object,
                                  cVector3d& a_segmentPointA,
                                  cVector3d& a_segmentPointB,
                                  cCollisionRecorder& a_recorder,
                                  cCollisionSettings& a_settings);

    //! This method renders the boxes of the tree at the display depth (see setDisplayDepth()).
    virtual void render(cRenderOptions& a_options);

    //! This method returns __true__ if a node is drawn at the current display depth.
    bool isNodeDisplayed(const cCollisionOBBNode& a_node) const;

    //! This method enables one color per tree level instead of m_color.
    void setColorByDepth(const bool a_colorByDepth) { m_colorByDepth = a_colorByDepth; }

    //! This method returns __true__ if the boxes are colored by tree level.
    bool getColorByDepth() const { return (m_colorByDepth); }

    //! This method sets the width of the lines of the boxes, in pixels.
    void setLineWidth(const double a_lineWidth) { m_lineWidth = a_lineWidth; }

    //! This method returns the width of the lines of the boxes, in pixels.
    double getLineWidth() const { return (m_lineWidth); }

    //! This method returns the depth of the deepest leaf of the tree.
    int getMaxDepth() const { return (m_tree.getStats().m_maxLeafDepth); }

    //! This method returns the OBB tree.
    const cCollisionOBBTree& getTree() const { return (m_tree); }

    //! This method returns the mesh of the detector.
    cMesh* getMesh() const { return (m_mesh); }


    //--------------------------------------------------------------------------
    // PROTECTED MEMBERS:
    //--------------------------------------------------------------------------

protected:

    //! If __true__, the boxes are colored by tree level.
    bool m_colorByDepth;

    //! Width of the lines of the boxes, in pixels.
    double m_lineWidth;

    //! Mesh of the detector.
    cMesh* m_mesh;

    //! OBB tree of the mesh (local coordinates of the mesh).
    cCollisionOBBTree m_tree;

    //! Settings used to build the tree.
    cOBBTreeSettings m_settings;
};


//------------------------------------------------------------------------------
/*!
    \brief
    Creates an OBB tree collision detector for a mesh and attaches it to the
    mesh (the previous detector of the mesh is deleted; the mesh owns the new
    one).

    \return The detector, or NULL if the mesh has no triangle.
*/
//------------------------------------------------------------------------------
cCollisionOBB* cCreateOBBCollisionDetector(cMesh* a_mesh,
                                           const cOBBTreeSettings& a_settings = cOBBTreeSettings());

//------------------------------------------------------------------------------
/*!
    \brief
    Tests whether two meshes intersect (e.g. a hand and an object).

    \details
    Both meshes must use a cCollisionOBB detector. The meshes are placed at
    their global poses (getGlobalRot(), getGlobalPos()), which CHAI3D updates
    in cWorld::computeGlobalPositions(). \n

    The returned pairs hold triangle indices of the meshes (index in
    cMesh::m_triangles).

    \param  a_meshA  First mesh.
    \param  a_meshB  Second mesh.
    \param  a_pairs  Output: all intersecting triangle pairs. If NULL, the
                     query stops at the first contact (fastest).
    \param  a_stats  Optional: counters incremented by the query.

    \return __true__ if the meshes intersect.
*/
//------------------------------------------------------------------------------
bool cComputeMeshIntersection(cMesh* a_meshA,
                              cMesh* a_meshB,
                              std::vector<cOBBTrianglePair>* a_pairs = NULL,
                              cOBBQueryStats* a_stats = NULL);

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
