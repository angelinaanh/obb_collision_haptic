//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBB.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CCollisionOBB.h"
//------------------------------------------------------------------------------
#include "CCollisionOBBDraw.h"
#include "system/CGlobals.h"
#include "world/CMesh.h"
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//==============================================================================
/*!
    Constructor of cCollisionOBB.
*/
//==============================================================================
cCollisionOBB::cCollisionOBB() : m_colorByDepth(false), m_lineWidth(1.0), m_mesh(NULL)
{
    // the tool radius is applied at query time (cCollisionSettings::m_collisionRadius)
    m_radiusAroundElements = 0.0;
}


//==============================================================================
/*!
    This method builds the OBB tree of a mesh.

    \param  a_mesh      Mesh.
    \param  a_settings  Construction settings.

    \return __true__ if the tree was built (the mesh has triangles).
*/
//==============================================================================
bool cCollisionOBB::initialize(cMesh* a_mesh, const cOBBTreeSettings& a_settings)
{
    m_mesh = a_mesh;
    m_settings = a_settings;
    return (m_tree.build(a_mesh, a_settings));
}


//==============================================================================
/*!
    This method rebuilds the tree from the current triangles of the mesh.
*/
//==============================================================================
void cCollisionOBB::update()
{
    m_tree.build(m_mesh, m_settings);
}


//==============================================================================
/*!
    This method computes the collisions between segment [A, B], expressed in
    the local frame of the mesh, and the triangles of the mesh.

    The tree returns the triangles of the leaves whose boxes, inflated by the
    collision radius, are crossed by the segment. Each candidate is then tested
    by the triangle array of the mesh, which records the collision events.
*/
//==============================================================================
bool cCollisionOBB::computeCollision(cGenericObject* a_object,
                                     cVector3d& a_segmentPointA,
                                     cVector3d& a_segmentPointB,
                                     cCollisionRecorder& a_recorder,
                                     cCollisionSettings& a_settings)
{
    if ((m_mesh == NULL) || m_tree.isEmpty())
    {
        return (false);
    }

    std::vector<unsigned int> candidates;
    m_tree.querySegment(a_segmentPointA, a_segmentPointB, a_settings.m_collisionRadius, candidates);

    bool hit = false;
    for (size_t i=0; i<candidates.size(); i++)
    {
        const int index = m_tree.getMeshTriangle(candidates[i]);
        if (m_mesh->m_triangles->computeCollision(index,
                                                  a_object,
                                                  a_segmentPointA,
                                                  a_segmentPointB,
                                                  a_recorder,
                                                  a_settings))
        {
            hit = true;
        }
    }

    return (hit);
}


//==============================================================================
/*!
    This method returns __true__ if a node is drawn at the current display
    depth d:
    - d < 0: nodes of depth 0 to |d|;
    - d >= 0: nodes of depth d, and leaves shallower than d.
*/
//==============================================================================
bool cCollisionOBB::isNodeDisplayed(const cCollisionOBBNode& a_node) const
{
    if (m_displayDepth < 0)
    {
        return (a_node.m_depth <= -m_displayDepth);
    }
    return ((a_node.m_depth == m_displayDepth) || (a_node.isLeaf() && (a_node.m_depth < m_displayDepth)));
}


//==============================================================================
/*!
    This method renders the boxes of the tree with OpenGL, as wireframes.

    CHAI3D calls this method while the local frame of the mesh is on the
    OpenGL matrix stack (see cGenericObject::renderSceneGraph()) and only when
    the display of the collision detector is enabled on the mesh
    (cGenericObject::setShowCollisionDetector()). The corners of the boxes are
    in the local frame of the mesh, so the boxes follow the mesh.
*/
//==============================================================================
void cCollisionOBB::render(cRenderOptions& /*a_options*/)
{
#ifdef C_USE_OPENGL

    const std::vector<cCollisionOBBNode>& nodes = m_tree.getNodes();
    if (nodes.empty())
    {
        return;
    }
    const double maxDepth = (double)cMax(1, getMaxDepth());

    glDisable(GL_LIGHTING);
    glLineWidth((GLfloat)m_lineWidth);
    glColor4fv(m_color.getData());

    glBegin(GL_LINES);
    for (size_t i=0; i<nodes.size(); i++)
    {
        if (!isNodeDisplayed(nodes[i])) continue;

        if (m_colorByDepth)
        {
            glColor4fv(cOBBDepthColor(nodes[i].m_depth / maxDepth).getData());
        }

        cVector3d corners[8];
        cGetOBBCorners(nodes[i].m_box, corners);
        for (int e=0; e<12; e++)
        {
            const cVector3d& a = corners[C_OBB_EDGES[e][0]];
            const cVector3d& b = corners[C_OBB_EDGES[e][1]];
            glVertex3d(a(0), a(1), a(2));
            glVertex3d(b(0), b(1), b(2));
        }
    }
    glEnd();

    glLineWidth(1.0f);
    glEnable(GL_LIGHTING);

#endif
}


//==============================================================================
/*!
    This function creates an OBB tree collision detector and attaches it to a
    mesh, replacing (and deleting) the previous detector of the mesh.
*/
//==============================================================================
cCollisionOBB* cCreateOBBCollisionDetector(cMesh* a_mesh, const cOBBTreeSettings& a_settings)
{
    if (a_mesh == NULL)
    {
        return (NULL);
    }

    cCollisionOBB* detector = new cCollisionOBB();
    if (!detector->initialize(a_mesh, a_settings))
    {
        delete detector;
        return (NULL);
    }

    delete a_mesh->getCollisionDetector();
    a_mesh->setCollisionDetector(detector);

    return (detector);
}


//==============================================================================
/*!
    This function tests whether two meshes with OBB detectors intersect at
    their global poses.
*/
//==============================================================================
bool cComputeMeshIntersection(cMesh* a_meshA,
                              cMesh* a_meshB,
                              std::vector<cOBBTrianglePair>* a_pairs,
                              cOBBQueryStats* a_stats)
{
    if (a_pairs != NULL) a_pairs->clear();
    if ((a_meshA == NULL) || (a_meshB == NULL)) return (false);

    const cCollisionOBB* detectorA = dynamic_cast<cCollisionOBB*>(a_meshA->getCollisionDetector());
    const cCollisionOBB* detectorB = dynamic_cast<cCollisionOBB*>(a_meshB->getCollisionDetector());
    if ((detectorA == NULL) || (detectorB == NULL)) return (false);

    const cCollisionOBBTree& treeA = detectorA->getTree();
    const cCollisionOBBTree& treeB = detectorB->getTree();

    const bool hit = cCollisionOBBTree::intersect(treeA, a_meshA->getGlobalRot(), a_meshA->getGlobalPos(),
                                                  treeB, a_meshB->getGlobalRot(), a_meshB->getGlobalPos(),
                                                  a_pairs, a_stats);

    // tree indices -> mesh triangle indices
    if (a_pairs != NULL)
    {
        for (size_t i=0; i<a_pairs->size(); i++)
        {
            (*a_pairs)[i].m_triangleA = (unsigned int)treeA.getMeshTriangle((*a_pairs)[i].m_triangleA);
            (*a_pairs)[i].m_triangleB = (unsigned int)treeB.getMeshTriangle((*a_pairs)[i].m_triangleB);
        }
    }

    return (hit);
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
