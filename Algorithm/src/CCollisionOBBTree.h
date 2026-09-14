//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBTree.h
    \ingroup    obb

    \brief
    Bounding volume hierarchy of oriented boxes built with the Surface Area
    Heuristic (SAH), with segment and tree-tree intersection queries.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef CCollisionOBBTreeH
#define CCollisionOBBTreeH
//------------------------------------------------------------------------------
#include "CCollisionOBBBox.h"
#include "CCollisionOBBRefine.h"
//------------------------------------------------------------------------------
#include <vector>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

class cMesh;

//! Strategy used to split a node in two.
enum cOBBSplitMethod
{
    //! Surface Area Heuristic: sweep of all split planes, minimum expected cost.
    C_OBB_SPLIT_SAH,

    //! Median of the triangle centroids along the longest box axis (reference).
    C_OBB_SPLIT_MEDIAN
};


//==============================================================================
/*!
    \struct     cOBBTreeSettings
    \ingroup    obb

    \brief
    Parameters of the OBB tree construction.

    \details
    __SAH__. For a node P with N triangles, every split plane orthogonal to
    one of the three axes of the node's OBB is evaluated (triangles sorted by
    centroid, all N-1 positions per axis). The cost of a split into L and R is

    \f[ Cost = C_{node} + \frac{S_L}{S_P} N_L + \frac{S_R}{S_P} N_R \f]

    where \f$S\f$ is a surface area and the cost of one triangle test is 1.
    \f$S_L/S_P\f$ is the probability that a segment crossing P also crosses L.
    The areas are those of the boxes aligned with the parent's OBB axes (they
    are updated incrementally during the sweep); the children then get their
    own tight OBB. The node becomes a leaf when no split is cheaper than
    testing its N triangles (\f$Cost \ge N\f$) and N is at most
    \ref m_maxLeafSize. \n

    __Balance__. Only splits leaving at least \ref m_minSplitFraction of the
    triangles on each side are considered, which bounds the depth of the tree
    (\f$O(\log N)\f$) while leaving SAH free to choose among the other planes. \n

    __Attention near the root__. Pruning a node near the root removes many
    triangles from consideration, so the boxes of the first
    \ref m_topLevels levels are refined with the wider and finer settings
    \ref m_topRefine; deeper boxes use \ref m_refine.
*/
//==============================================================================
struct cOBBTreeSettings
{
    //! Constructor of cOBBTreeSettings.
    cOBBTreeSettings() :
        m_splitMethod(C_OBB_SPLIT_SAH),
        m_nodeCost(1.0),
        m_maxLeafSize(8),
        m_minSplitFraction(0.1),
        m_maxDepth(64),
        m_topLevels(4)
    {
        m_topRefine.m_maxIterations = 10;
        m_topRefine.m_anglesDeg.clear();
        m_topRefine.m_anglesDeg.push_back(15.0);
        m_topRefine.m_anglesDeg.push_back(5.0);
        m_topRefine.m_anglesDeg.push_back(1.0);
        m_topRefine.m_anglesDeg.push_back(0.2);
    }

    //! Split strategy.
    cOBBSplitMethod m_splitMethod;

    //! Cost C_node of testing a box, relative to the cost of testing a triangle.
    double m_nodeCost;

    //! SAH may stop splitting nodes with at most this many triangles (nodes with more are always split).
    int m_maxLeafSize;

    //! Minimum fraction of the triangles of a node that each child receives (SAH).
    double m_minSplitFraction;

    //! Maximum depth of the tree.
    int m_maxDepth;

    //! Number of levels below the root that use m_topRefine.
    int m_topLevels;

    //! Refinement of the boxes of the first m_topLevels levels.
    cOBBRefineSettings m_topRefine;

    //! Refinement of the other boxes.
    cOBBRefineSettings m_refine;
};


//==============================================================================
/*!
    \struct     cCollisionOBBNode
    \ingroup    obb

    \brief
    Node of an OBB tree.

    \details
    Internal nodes have two children. Leaves reference a range of triangles
    in the leaf order of the tree (see cCollisionOBBTree::getLeafTriangle()).
*/
//==============================================================================
struct cCollisionOBBNode
{
    //! Constructor of cCollisionOBBNode.
    cCollisionOBBNode() : m_left(-1), m_right(-1), m_first(0), m_count(0), m_depth(0) {}

    //! __true__ if the node has no children.
    inline bool isLeaf() const { return (m_left < 0); }

    //! Tight oriented box around the triangles of the node.
    cCollisionOBBBox m_box;

    //! Index of the left child (-1 for a leaf).
    int m_left;

    //! Index of the right child (-1 for a leaf).
    int m_right;

    //! First triangle of the node in the leaf order.
    unsigned int m_first;

    //! Number of triangles below the node.
    unsigned int m_count;

    //! Depth of the node (root = 0).
    int m_depth;
};


//==============================================================================
/*!
    \struct     cOBBTreeStats
    \ingroup    obb

    \brief
    Statistics of an OBB tree, computed after construction.
*/
//==============================================================================
struct cOBBTreeStats
{
    cOBBTreeStats() :
        m_numTriangles(0), m_numNodes(0), m_numLeaves(0),
        m_minLeafDepth(0), m_maxLeafDepth(0), m_meanLeafDepth(0.0),
        m_maxLeafSize(0), m_meanLeafSize(0.0),
        m_sahCost(0.0), m_volumeRatio(0.0), m_meanSiblingOverlap(0.0),
        m_buildTime(0.0) {}

    int m_numTriangles;
    int m_numNodes;
    int m_numLeaves;
    int m_minLeafDepth;
    int m_maxLeafDepth;
    double m_meanLeafDepth;
    int m_maxLeafSize;
    double m_meanLeafSize;

    //! Expected cost of a segment query crossing the root, from the actual OBB areas.
    double m_sahCost;

    //! Sum of the volumes of all nodes divided by the volume of the root.
    double m_volumeRatio;

    //! Mean overlap volume of sibling boxes (in the parent frame) relative to the parent box.
    double m_meanSiblingOverlap;

    //! Construction time in seconds.
    double m_buildTime;
};


//! Nearest intersection between a segment and the triangles of a tree.
struct cOBBSegmentHit
{
    cOBBSegmentHit() : m_hit(false), m_t(0.0), m_triangle(0), m_meshTriangle(-1) {}

    //! __true__ if the segment intersects a triangle.
    bool m_hit;

    //! Segment parameter in [0, 1] of the intersection point.
    double m_t;

    //! Intersection point.
    cVector3d m_point;

    //! Unit normal of the intersected triangle.
    cVector3d m_normal;

    //! Index of the triangle in the tree.
    unsigned int m_triangle;

    //! Index of the triangle in the source mesh (-1 if not built from a mesh).
    int m_meshTriangle;
};


//! Pair of intersecting triangles (indices in tree A and tree B).
struct cOBBTrianglePair
{
    unsigned int m_triangleA;
    unsigned int m_triangleB;
};


//! Cost counters of a query.
struct cOBBQueryStats
{
    cOBBQueryStats() : m_boxTests(0), m_triangleTests(0) {}

    //! Number of box tests (segment-box or box-box).
    long long m_boxTests;

    //! Number of triangle tests (segment-triangle or triangle-triangle).
    long long m_triangleTests;
};


//==============================================================================
/*!
    \class      cCollisionOBBTree
    \ingroup    obb

    \brief
    Bounding volume hierarchy of tight oriented bounding boxes.

    \details
    Every node holds the tight OBB of its triangles (convex hull, continuous
    PCA and iterative refinement, see cComputeOBB()). The hierarchy is built
    top-down with the Surface Area Heuristic (see cOBBTreeSettings). \n

    All coordinates are expressed in the local frame of the geometry (the
    local frame of the mesh when built from a cMesh).
*/
//==============================================================================
class cCollisionOBBTree
{
    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of cCollisionOBBTree.
    cCollisionOBBTree() : m_overlapSum(0.0), m_numSplits(0) {}

    //! Destructor of cCollisionOBBTree.
    virtual ~cCollisionOBBTree() {}


    //--------------------------------------------------------------------------
    // PUBLIC METHODS - CONSTRUCTION:
    //--------------------------------------------------------------------------

public:

    //! This method builds the tree of a triangle soup (3 vertex indices per triangle).
    bool build(const std::vector<cVector3d>& a_vertices,
               const std::vector<unsigned int>& a_triangles,
               const cOBBTreeSettings& a_settings = cOBBTreeSettings());

    //! This method builds the tree of the allocated triangles of a mesh, in the local frame of the mesh.
    bool build(cMesh* a_mesh,
               const cOBBTreeSettings& a_settings = cOBBTreeSettings());

    //! This method clears the tree.
    void clear();


    //--------------------------------------------------------------------------
    // PUBLIC METHODS - ACCESS:
    //--------------------------------------------------------------------------

public:

    //! This method returns __true__ if the tree has no triangle.
    bool isEmpty() const { return (m_nodes.empty()); }

    //! This method returns the nodes of the tree (root = 0).
    const std::vector<cCollisionOBBNode>& getNodes() const { return (m_nodes); }

    //! This method returns the number of triangles.
    unsigned int getNumTriangles() const { return ((unsigned int)m_triangles.size() / 3); }

    //! This method returns the triangle stored at a position of the leaf order.
    unsigned int getLeafTriangle(unsigned int a_position) const { return (m_order[a_position]); }

    //! This method returns the vertices of a triangle.
    void getTriangle(unsigned int a_triangle, cVector3d& a_v0, cVector3d& a_v1, cVector3d& a_v2) const
    {
        a_v0 = m_vertices[m_triangles[3*a_triangle+0]];
        a_v1 = m_vertices[m_triangles[3*a_triangle+1]];
        a_v2 = m_vertices[m_triangles[3*a_triangle+2]];
    }

    //! This method returns the index of a triangle in the source mesh (-1 if not built from a mesh).
    int getMeshTriangle(unsigned int a_triangle) const
    {
        return (m_meshTriangles.empty() ? -1 : m_meshTriangles[a_triangle]);
    }

    //! This method returns the settings used for the last construction.
    const cOBBTreeSettings& getSettings() const { return (m_settings); }

    //! This method returns the statistics of the tree.
    const cOBBTreeStats& getStats() const { return (m_stats); }


    //--------------------------------------------------------------------------
    // PUBLIC METHODS - QUERIES:
    //--------------------------------------------------------------------------

public:

    //! This method finds the nearest intersection between segment [A, B] and the triangles.
    bool intersectSegment(const cVector3d& a_pointA,
                          const cVector3d& a_pointB,
                          cOBBSegmentHit& a_hit,
                          cOBBQueryStats* a_stats = NULL) const;

    //! This method collects the triangles of all leaves whose (inflated) boxes are crossed by segment [A, B].
    void querySegment(const cVector3d& a_pointA,
                      const cVector3d& a_pointB,
                      double a_radius,
                      std::vector<unsigned int>& a_triangles,
                      cOBBQueryStats* a_stats = NULL) const;

    //! This method finds the intersecting triangles of two trees placed at given poses.
    static bool intersect(const cCollisionOBBTree& a_treeA,
                          const cMatrix3d& a_rotA,
                          const cVector3d& a_posA,
                          const cCollisionOBBTree& a_treeB,
                          const cMatrix3d& a_rotB,
                          const cVector3d& a_posB,
                          std::vector<cOBBTrianglePair>* a_pairs = NULL,
                          cOBBQueryStats* a_stats = NULL);


    //--------------------------------------------------------------------------
    // PROTECTED METHODS:
    //--------------------------------------------------------------------------

protected:

    //! This method builds the node for triangles [a_first, a_first + a_count) of the leaf order.
    int buildNode(unsigned int a_first, unsigned int a_count, int a_depth);

    //! This method reorders the triangles of a node with SAH; returns the size of the left child (0: leaf).
    unsigned int splitSAH(unsigned int a_first, unsigned int a_count, const cCollisionOBBBox& a_box);

    //! This method reorders the triangles of a node at the median; returns the size of the left child.
    unsigned int splitMedian(unsigned int a_first, unsigned int a_count, const cCollisionOBBBox& a_box);

    //! This method computes the statistics of the tree.
    void computeStats();


    //--------------------------------------------------------------------------
    // PROTECTED MEMBERS:
    //--------------------------------------------------------------------------

protected:

    //! Vertices (local coordinates).
    std::vector<cVector3d> m_vertices;

    //! Triangles (3 vertex indices per triangle).
    std::vector<unsigned int> m_triangles;

    //! Index of each triangle in the source mesh.
    std::vector<int> m_meshTriangles;

    //! Triangle indices in leaf order: the triangles of every node are contiguous.
    std::vector<unsigned int> m_order;

    //! Nodes (root = 0).
    std::vector<cCollisionOBBNode> m_nodes;

    //! Settings of the last construction.
    cOBBTreeSettings m_settings;

    //! Statistics of the tree.
    cOBBTreeStats m_stats;

    //! Construction scratch: last node that used each vertex.
    std::vector<unsigned int> m_vertexStamp;

    //! Construction scratch: sum of the sibling overlap ratios.
    double m_overlapSum;

    //! Construction scratch: number of splits.
    int m_numSplits;
};

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
