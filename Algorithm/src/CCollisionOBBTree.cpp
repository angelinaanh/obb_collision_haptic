//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBTree.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CCollisionOBBTree.h"
//------------------------------------------------------------------------------
#include "CCollisionOBBFit.h"
#include "CCollisionOBBIntersect.h"
#include "timers/CPrecisionClock.h"
#include "world/CMesh.h"
//------------------------------------------------------------------------------
#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <utility>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace {
//------------------------------------------------------------------------------

//! Axis-aligned bounds in the frame of a node (the axes of its OBB).
struct cFrameBounds
{
    double m_lo[3];
    double m_hi[3];

    cFrameBounds()
    {
        for (int k=0; k<3; k++) { m_lo[k] = DBL_MAX; m_hi[k] = -DBL_MAX; }
    }

    void enclose(const double* a_lo, const double* a_hi)
    {
        for (int k=0; k<3; k++)
        {
            if (a_lo[k] < m_lo[k]) m_lo[k] = a_lo[k];
            if (a_hi[k] > m_hi[k]) m_hi[k] = a_hi[k];
        }
    }

    double size(int a_axis) const { return (cMax(0.0, m_hi[a_axis] - m_lo[a_axis])); }

    double area() const
    {
        const double x = size(0), y = size(1), z = size(2);
        return (2.0 * (x * y + y * z + z * x));
    }

    double volume() const { return (size(0) * size(1) * size(2)); }
};


//! Volume of the intersection of two bounds.
double overlapVolume(const cFrameBounds& a_a, const cFrameBounds& a_b)
{
    double volume = 1.0;
    for (int k=0; k<3; k++)
    {
        const double o = cMin(a_a.m_hi[k], a_b.m_hi[k]) - cMax(a_a.m_lo[k], a_b.m_lo[k]);
        if (o <= 0.0) return (0.0);
        volume *= o;
    }
    return (volume);
}


//! Projections of the triangles of a node on the axes of its box.
struct cTriangleProjections
{
    std::vector<double> m_lo;       // 3 per triangle: min along each axis
    std::vector<double> m_hi;       // 3 per triangle: max along each axis
    std::vector<double> m_centroid; // 3 per triangle: centroid along each axis

    const double* lo(unsigned int a_i) const { return (&m_lo[3*a_i]); }
    const double* hi(unsigned int a_i) const { return (&m_hi[3*a_i]); }
    double centroid(unsigned int a_i, int a_axis) const { return (m_centroid[3*a_i + a_axis]); }
};


//! Sorts triangle slots 0..n-1 by centroid along an axis (ties by slot, for determinism).
void sortByCentroid(const cTriangleProjections& a_proj, int a_axis, std::vector<unsigned int>& a_order)
{
    for (unsigned int i=0; i<(unsigned int)a_order.size(); i++) a_order[i] = i;
    std::sort(a_order.begin(), a_order.end(), [&](unsigned int a, unsigned int b)
    {
        const double ca = a_proj.centroid(a, a_axis);
        const double cb = a_proj.centroid(b, a_axis);
        return ((ca < cb) || ((ca == cb) && (a < b)));
    });
}


//! Volume of a stored box.
inline double boxVolume(const cCollisionOBBBox& a_box)
{
    return (8.0 * (double)a_box.m_halfExtents[0] * a_box.m_halfExtents[1] * a_box.m_halfExtents[2]);
}


//! Surface area of a stored box.
inline double boxArea(const cCollisionOBBBox& a_box)
{
    const double x = a_box.m_halfExtents[0], y = a_box.m_halfExtents[1], z = a_box.m_halfExtents[2];
    return (8.0 * (x * y + y * z + z * x));
}

//------------------------------------------------------------------------------
} // anonymous namespace
//------------------------------------------------------------------------------


//==============================================================================
/*!
    This method clears the tree.
*/
//==============================================================================
void cCollisionOBBTree::clear()
{
    m_vertices.clear();
    m_triangles.clear();
    m_meshTriangles.clear();
    m_order.clear();
    m_nodes.clear();
    m_stats = cOBBTreeStats();
    m_overlapSum = 0.0;
    m_numSplits = 0;
}


//==============================================================================
/*!
    This method builds the tree of a triangle soup.

    \param  a_vertices   Vertices.
    \param  a_triangles  Triangles (3 vertex indices per triangle).
    \param  a_settings   Construction settings.

    \return __false__ if there is no triangle or an index is out of range.
*/
//==============================================================================
bool cCollisionOBBTree::build(const std::vector<cVector3d>& a_vertices,
                              const std::vector<unsigned int>& a_triangles,
                              const cOBBTreeSettings& a_settings)
{
    cPrecisionClock clock;
    clock.start(true);

    clear();
    m_settings = a_settings;

    const unsigned int numTriangles = (unsigned int)(a_triangles.size() / 3);
    for (size_t i=0; i<3*(size_t)numTriangles; i++)
    {
        if (a_triangles[i] >= a_vertices.size()) return (false);
    }
    if (numTriangles == 0) return (false);

    m_vertices = a_vertices;
    m_triangles.assign(a_triangles.begin(), a_triangles.begin() + 3 * (size_t)numTriangles);
    m_order.resize(numTriangles);
    for (unsigned int i=0; i<numTriangles; i++) m_order[i] = i;

    m_vertexStamp.assign(m_vertices.size(), 0);
    m_nodes.reserve(2 * (size_t)numTriangles);

    buildNode(0, numTriangles, 0);

    std::vector<unsigned int>().swap(m_vertexStamp);
    computeStats();
    m_stats.m_buildTime = clock.getCurrentTimeSeconds();

    return (true);
}


//==============================================================================
/*!
    This method builds the tree of the allocated triangles of a mesh. The tree
    is expressed in the local frame of the mesh.
*/
//==============================================================================
bool cCollisionOBBTree::build(cMesh* a_mesh, const cOBBTreeSettings& a_settings)
{
    if (a_mesh == NULL)
    {
        clear();
        return (false);
    }

    cVertexArrayPtr vertexArray = a_mesh->m_vertices;
    cTriangleArrayPtr triangleArray = a_mesh->m_triangles;

    std::vector<cVector3d> vertices(vertexArray->getNumElements());
    for (unsigned int i=0; i<(unsigned int)vertices.size(); i++)
    {
        vertices[i] = vertexArray->getLocalPos(i);
    }

    std::vector<unsigned int> triangles;
    std::vector<int> meshTriangles;
    const unsigned int numTriangles = triangleArray->getNumElements();
    for (unsigned int i=0; i<numTriangles; i++)
    {
        if (!triangleArray->getAllocated(i)) continue;
        triangles.push_back(triangleArray->getVertexIndex0(i));
        triangles.push_back(triangleArray->getVertexIndex1(i));
        triangles.push_back(triangleArray->getVertexIndex2(i));
        meshTriangles.push_back((int)i);
    }

    if (!build(vertices, triangles, a_settings))
    {
        return (false);
    }
    m_meshTriangles.swap(meshTriangles);
    return (true);
}


//==============================================================================
/*!
    This method builds the node that holds triangles [a_first, a_first + a_count)
    of the leaf order, then its children.

    \return Index of the node.
*/
//==============================================================================
int cCollisionOBBTree::buildNode(unsigned int a_first, unsigned int a_count, int a_depth)
{
    const int index = (int)m_nodes.size();
    m_nodes.push_back(cCollisionOBBNode());
    m_nodes[index].m_first = a_first;
    m_nodes[index].m_count = a_count;
    m_nodes[index].m_depth = a_depth;

    // tight OBB of the vertices used by the triangles of the node
    cCollisionOBBBox box;
    {
        const unsigned int stamp = (unsigned int)index + 1;
        std::vector<cVector3d> points;
        for (unsigned int i=0; i<a_count; i++)
        {
            const unsigned int tri = m_order[a_first + i];
            for (int k=0; k<3; k++)
            {
                const unsigned int v = m_triangles[3*tri + k];
                if (m_vertexStamp[v] != stamp)
                {
                    m_vertexStamp[v] = stamp;
                    points.push_back(m_vertices[v]);
                }
            }
        }

        const cOBBRefineSettings& refine = (a_depth < m_settings.m_topLevels) ?
                                           m_settings.m_topRefine : m_settings.m_refine;
        box = cComputeOBB(points, NULL, refine);
    }
    m_nodes[index].m_box = box;

    if ((a_count <= 1) || (a_depth >= m_settings.m_maxDepth))
    {
        return (index);
    }

    const unsigned int split = (m_settings.m_splitMethod == C_OBB_SPLIT_SAH) ?
                               splitSAH(a_first, a_count, box) :
                               splitMedian(a_first, a_count, box);
    if (split == 0)
    {
        return (index);
    }

    // m_nodes may be reallocated by the recursion: store indices only
    const int left = buildNode(a_first, split, a_depth + 1);
    const int right = buildNode(a_first + split, a_count - split, a_depth + 1);
    m_nodes[index].m_left = left;
    m_nodes[index].m_right = right;

    return (index);
}


//------------------------------------------------------------------------------
// Projects the triangles [a_first, a_first + a_count) of the leaf order on the
// axes of a box.
//------------------------------------------------------------------------------
static void projectTriangles(const std::vector<cVector3d>& a_vertices,
                             const std::vector<unsigned int>& a_triangles,
                             const std::vector<unsigned int>& a_order,
                             unsigned int a_first,
                             unsigned int a_count,
                             const cCollisionOBBBox& a_box,
                             cTriangleProjections& a_proj)
{
    const cMatrix3d rot = a_box.getRotationMatrix();
    const cVector3d axis[3] = { rot.getCol0(), rot.getCol1(), rot.getCol2() };

    a_proj.m_lo.resize(3 * (size_t)a_count);
    a_proj.m_hi.resize(3 * (size_t)a_count);
    a_proj.m_centroid.resize(3 * (size_t)a_count);

    for (unsigned int i=0; i<a_count; i++)
    {
        const unsigned int tri = a_order[a_first + i];
        const cVector3d& v0 = a_vertices[a_triangles[3*tri + 0]];
        const cVector3d& v1 = a_vertices[a_triangles[3*tri + 1]];
        const cVector3d& v2 = a_vertices[a_triangles[3*tri + 2]];
        for (int k=0; k<3; k++)
        {
            const double t0 = cDot(v0, axis[k]);
            const double t1 = cDot(v1, axis[k]);
            const double t2 = cDot(v2, axis[k]);
            a_proj.m_lo[3*i + k] = cMin(t0, cMin(t1, t2));
            a_proj.m_hi[3*i + k] = cMax(t0, cMax(t1, t2));
            a_proj.m_centroid[3*i + k] = (t0 + t1 + t2) / 3.0;
        }
    }
}


//------------------------------------------------------------------------------
// Applies a split: a_order lists the slots of the node's triangles, the first
// a_split go to the left child. Also returns the sibling overlap ratio.
//------------------------------------------------------------------------------
static double applySplit(std::vector<unsigned int>& a_leafOrder,
                         unsigned int a_first,
                         const std::vector<unsigned int>& a_order,
                         unsigned int a_split,
                         const cTriangleProjections& a_proj)
{
    const unsigned int n = (unsigned int)a_order.size();

    cFrameBounds left, right, parent;
    for (unsigned int i=0; i<n; i++)
    {
        const unsigned int s = a_order[i];
        if (i < a_split) left.enclose(a_proj.lo(s), a_proj.hi(s));
        else             right.enclose(a_proj.lo(s), a_proj.hi(s));
        parent.enclose(a_proj.lo(s), a_proj.hi(s));
    }

    std::vector<unsigned int> reordered(n);
    for (unsigned int i=0; i<n; i++) reordered[i] = a_leafOrder[a_first + a_order[i]];
    std::copy(reordered.begin(), reordered.end(), a_leafOrder.begin() + a_first);

    const double volume = parent.volume();
    return ((volume > 0.0) ? overlapVolume(left, right) / volume : 0.0);
}


//==============================================================================
/*!
    This method evaluates the SAH cost of every split plane orthogonal to the
    axes of the node's box and applies the cheapest one.

    \return Number of triangles of the left child, or 0 if the node must stay
            a leaf.
*/
//==============================================================================
unsigned int cCollisionOBBTree::splitSAH(unsigned int a_first,
                                         unsigned int a_count,
                                         const cCollisionOBBBox& a_box)
{
    const unsigned int n = a_count;

    cTriangleProjections proj;
    projectTriangles(m_vertices, m_triangles, m_order, a_first, a_count, a_box, proj);

    // area of the parent, measured with the same kind of box as the children
    cFrameBounds parent;
    for (unsigned int i=0; i<n; i++) parent.enclose(proj.lo(i), proj.hi(i));
    const double parentArea = parent.area();

    // balance: each child receives at least this many triangles
    unsigned int minCount = (unsigned int)ceil(m_settings.m_minSplitFraction * n);
    minCount = cClamp(minCount, 1u, n / 2);

    std::vector<unsigned int> order(n), bestOrder;
    std::vector<double> rightArea(n, 0.0);
    double bestCost = DBL_MAX;
    unsigned int bestSplit = 0;

    for (int axis=0; axis<3; axis++)
    {
        sortByCentroid(proj, axis, order);

        // sweep from the right: area of the triangles order[j..n-1]
        cFrameBounds right;
        for (unsigned int j=n-1; j>=1; j--)
        {
            right.enclose(proj.lo(order[j]), proj.hi(order[j]));
            rightArea[j] = right.area();
        }

        // sweep from the left: j triangles on the left side
        bool improved = false;
        cFrameBounds left;
        for (unsigned int j=1; j<n; j++)
        {
            left.enclose(proj.lo(order[j-1]), proj.hi(order[j-1]));
            if ((j < minCount) || (n - j < minCount)) continue;

            // Cost = C_node + (S_L / S_P) N_L + (S_R / S_P) N_R
            double cost = m_settings.m_nodeCost;
            if (parentArea > 0.0)
            {
                cost += (left.area() * j + rightArea[j] * (n - j)) / parentArea;
            }
            else
            {
                cost += n;  // degenerate node: no area information
            }

            // on equal cost, prefer the more balanced split
            const bool better = (cost < bestCost) ||
                                ((cost == bestCost) && (abs((int)(2 * j) - (int)n) < abs((int)(2 * bestSplit) - (int)n)));
            if (better)
            {
                bestCost = cost;
                bestSplit = j;
                improved = true;
            }
        }

        if (improved)
        {
            bestOrder = order;
        }
    }

    // keep a leaf when no split is cheaper than testing all its triangles
    if ((bestSplit == 0) || ((n <= (unsigned int)m_settings.m_maxLeafSize) && (bestCost >= (double)n)))
    {
        return (0);
    }

    m_overlapSum += applySplit(m_order, a_first, bestOrder, bestSplit, proj);
    m_numSplits++;

    return (bestSplit);
}


//==============================================================================
/*!
    This method splits the triangles of a node at the median of their
    centroids along the longest axis of the node's box (reference strategy,
    down to one triangle per leaf).

    \return Number of triangles of the left child.
*/
//==============================================================================
unsigned int cCollisionOBBTree::splitMedian(unsigned int a_first,
                                            unsigned int a_count,
                                            const cCollisionOBBBox& a_box)
{
    cTriangleProjections proj;
    projectTriangles(m_vertices, m_triangles, m_order, a_first, a_count, a_box, proj);

    int axis = 0;
    for (int k=1; k<3; k++)
    {
        if (a_box.m_halfExtents[k] > a_box.m_halfExtents[axis]) axis = k;
    }

    std::vector<unsigned int> order(a_count);
    sortByCentroid(proj, axis, order);

    const unsigned int split = a_count / 2;
    m_overlapSum += applySplit(m_order, a_first, order, split, proj);
    m_numSplits++;

    return (split);
}


//==============================================================================
/*!
    This method computes the statistics of the tree.
*/
//==============================================================================
void cCollisionOBBTree::computeStats()
{
    m_stats = cOBBTreeStats();
    if (m_nodes.empty()) return;

    m_stats.m_numTriangles = (int)getNumTriangles();
    m_stats.m_numNodes = (int)m_nodes.size();
    m_stats.m_minLeafDepth = INT_MAX;

    double sah = 0.0;
    double volumeSum = 0.0;
    double depthSum = 0.0;
    for (size_t i=0; i<m_nodes.size(); i++)
    {
        const cCollisionOBBNode& node = m_nodes[i];
        const double area = boxArea(node.m_box);
        volumeSum += boxVolume(node.m_box);

        if (node.isLeaf())
        {
            m_stats.m_numLeaves++;
            m_stats.m_minLeafDepth = cMin(m_stats.m_minLeafDepth, node.m_depth);
            m_stats.m_maxLeafDepth = cMax(m_stats.m_maxLeafDepth, node.m_depth);
            m_stats.m_maxLeafSize = cMax(m_stats.m_maxLeafSize, (int)node.m_count);
            depthSum += node.m_depth;
            sah += area * node.m_count;
        }
        else
        {
            sah += area * m_settings.m_nodeCost;
        }
    }

    const double rootArea = boxArea(m_nodes[0].m_box);
    const double rootVolume = boxVolume(m_nodes[0].m_box);
    m_stats.m_meanLeafDepth = depthSum / m_stats.m_numLeaves;
    m_stats.m_meanLeafSize = (double)m_stats.m_numTriangles / m_stats.m_numLeaves;
    m_stats.m_sahCost = (rootArea > 0.0) ? sah / rootArea : 0.0;
    m_stats.m_volumeRatio = (rootVolume > 0.0) ? volumeSum / rootVolume : 0.0;
    m_stats.m_meanSiblingOverlap = (m_numSplits > 0) ? m_overlapSum / m_numSplits : 0.0;
}


//==============================================================================
/*!
    This method finds the nearest intersection between segment [A, B] and the
    triangles of the tree. Children are visited front to back, and a node is
    skipped when the segment enters its box beyond the nearest hit found so far.

    \param  a_pointA  Start point of the segment.
    \param  a_pointB  End point of the segment.
    \param  a_hit     Output: nearest intersection.
    \param  a_stats   Optional: counters incremented by the query.

    \return __true__ if the segment intersects a triangle.
*/
//==============================================================================
bool cCollisionOBBTree::intersectSegment(const cVector3d& a_pointA,
                                         const cVector3d& a_pointB,
                                         cOBBSegmentHit& a_hit,
                                         cOBBQueryStats* a_stats) const
{
    a_hit = cOBBSegmentHit();
    if (m_nodes.empty()) return (false);

    long long boxTests = 1;
    long long triangleTests = 0;
    double bestT = DBL_MAX;
    unsigned int bestTriangle = 0;

    std::vector<std::pair<int, double> > stack;
    double tRoot;
    if (cTestSegmentOBB(a_pointA, a_pointB, cOBBShape(m_nodes[0].m_box), 0.0, tRoot))
    {
        stack.push_back(std::make_pair(0, tRoot));
    }

    while (!stack.empty())
    {
        const std::pair<int, double> entry = stack.back();
        stack.pop_back();
        if (entry.second > bestT) continue;

        const cCollisionOBBNode& node = m_nodes[entry.first];
        if (node.isLeaf())
        {
            for (unsigned int i=0; i<node.m_count; i++)
            {
                const unsigned int tri = m_order[node.m_first + i];
                cVector3d v0, v1, v2;
                getTriangle(tri, v0, v1, v2);
                double t;
                triangleTests++;
                if (cTestSegmentTriangle(a_pointA, a_pointB, v0, v1, v2, t) && (t < bestT))
                {
                    bestT = t;
                    bestTriangle = tri;
                }
            }
            continue;
        }

        double tLeft, tRight;
        const bool hitLeft = cTestSegmentOBB(a_pointA, a_pointB, cOBBShape(m_nodes[node.m_left].m_box), 0.0, tLeft);
        const bool hitRight = cTestSegmentOBB(a_pointA, a_pointB, cOBBShape(m_nodes[node.m_right].m_box), 0.0, tRight);
        boxTests += 2;

        // push the farther child first so that the nearer one is visited first
        if (hitLeft && hitRight)
        {
            if (tLeft <= tRight)
            {
                stack.push_back(std::make_pair(node.m_right, tRight));
                stack.push_back(std::make_pair(node.m_left, tLeft));
            }
            else
            {
                stack.push_back(std::make_pair(node.m_left, tLeft));
                stack.push_back(std::make_pair(node.m_right, tRight));
            }
        }
        else if (hitLeft)
        {
            stack.push_back(std::make_pair(node.m_left, tLeft));
        }
        else if (hitRight)
        {
            stack.push_back(std::make_pair(node.m_right, tRight));
        }
    }

    if (a_stats != NULL)
    {
        a_stats->m_boxTests += boxTests;
        a_stats->m_triangleTests += triangleTests;
    }

    if (bestT == DBL_MAX)
    {
        return (false);
    }

    cVector3d v0, v1, v2;
    getTriangle(bestTriangle, v0, v1, v2);
    a_hit.m_hit = true;
    a_hit.m_t = bestT;
    a_hit.m_point = a_pointA + bestT * (a_pointB - a_pointA);
    a_hit.m_normal = cNormalize(cCross(v1 - v0, v2 - v0));
    a_hit.m_triangle = bestTriangle;
    a_hit.m_meshTriangle = getMeshTriangle(bestTriangle);

    return (true);
}


//==============================================================================
/*!
    This method collects the triangles of all leaves whose boxes, inflated by
    a_radius, are crossed by segment [A, B]. It is the broad phase used by the
    CHAI3D collision detector (cCollisionOBB).
*/
//==============================================================================
void cCollisionOBBTree::querySegment(const cVector3d& a_pointA,
                                     const cVector3d& a_pointB,
                                     double a_radius,
                                     std::vector<unsigned int>& a_triangles,
                                     cOBBQueryStats* a_stats) const
{
    a_triangles.clear();
    if (m_nodes.empty()) return;

    long long boxTests = 0;
    std::vector<int> stack(1, 0);
    while (!stack.empty())
    {
        const int index = stack.back();
        stack.pop_back();

        const cCollisionOBBNode& node = m_nodes[index];
        double t;
        boxTests++;
        if (!cTestSegmentOBB(a_pointA, a_pointB, cOBBShape(node.m_box), a_radius, t)) continue;

        if (node.isLeaf())
        {
            for (unsigned int i=0; i<node.m_count; i++)
            {
                a_triangles.push_back(m_order[node.m_first + i]);
            }
        }
        else
        {
            stack.push_back(node.m_right);
            stack.push_back(node.m_left);
        }
    }

    if (a_stats != NULL)
    {
        a_stats->m_boxTests += boxTests;
    }
}


//==============================================================================
/*!
    This method finds the pairs of intersecting triangles of two trees placed
    at given poses (e.g. a hand and an object). Pairs of nodes are traversed
    simultaneously; a pair is pruned as soon as its boxes are separated, and
    the node with the larger box is descended first.

    \param  a_treeA  First tree.
    \param  a_rotA   Rotation of the first tree (local to world).
    \param  a_posA   Position of the first tree (local to world).
    \param  a_treeB  Second tree.
    \param  a_rotB   Rotation of the second tree (local to world).
    \param  a_posB   Position of the second tree (local to world).
    \param  a_pairs  Output: all intersecting triangle pairs. If NULL, the query
                     stops at the first intersection.
    \param  a_stats  Optional: counters incremented by the query.

    \return __true__ if the trees intersect.
*/
//==============================================================================
bool cCollisionOBBTree::intersect(const cCollisionOBBTree& a_treeA,
                                  const cMatrix3d& a_rotA,
                                  const cVector3d& a_posA,
                                  const cCollisionOBBTree& a_treeB,
                                  const cMatrix3d& a_rotB,
                                  const cVector3d& a_posB,
                                  std::vector<cOBBTrianglePair>* a_pairs,
                                  cOBBQueryStats* a_stats)
{
    if (a_pairs != NULL) a_pairs->clear();
    if (a_treeA.isEmpty() || a_treeB.isEmpty()) return (false);

    // pose of tree B in the frame of tree A
    const cMatrix3d rotAT = cTranspose(a_rotA);
    const cMatrix3d R = rotAT * a_rotB;
    const cVector3d t = rotAT * (a_posB - a_posA);

    long long boxTests = 0;
    long long triangleTests = 0;
    bool found = false;

    std::vector<std::pair<int, int> > stack(1, std::make_pair(0, 0));
    std::vector<cVector3d> trianglesB;
    std::vector<unsigned int> indicesB;

    while (!stack.empty())
    {
        const std::pair<int, int> pair = stack.back();
        stack.pop_back();

        const cCollisionOBBNode& nodeA = a_treeA.m_nodes[pair.first];
        const cCollisionOBBNode& nodeB = a_treeB.m_nodes[pair.second];

        boxTests++;
        if (!cTestOBBOBB(cOBBShape(nodeA.m_box), cOBBShape(nodeB.m_box).transformed(R, t)))
        {
            continue;
        }

        if (nodeA.isLeaf() && nodeB.isLeaf())
        {
            // triangles of B in the frame of A
            trianglesB.resize(3 * (size_t)nodeB.m_count);
            indicesB.resize(nodeB.m_count);
            for (unsigned int j=0; j<nodeB.m_count; j++)
            {
                indicesB[j] = a_treeB.m_order[nodeB.m_first + j];
                cVector3d q0, q1, q2;
                a_treeB.getTriangle(indicesB[j], q0, q1, q2);
                trianglesB[3*j + 0] = R * q0 + t;
                trianglesB[3*j + 1] = R * q1 + t;
                trianglesB[3*j + 2] = R * q2 + t;
            }

            for (unsigned int i=0; i<nodeA.m_count; i++)
            {
                const unsigned int triA = a_treeA.m_order[nodeA.m_first + i];
                cVector3d p0, p1, p2;
                a_treeA.getTriangle(triA, p0, p1, p2);

                for (unsigned int j=0; j<nodeB.m_count; j++)
                {
                    triangleTests++;
                    if (cTestTriangleTriangle(p0, p1, p2, trianglesB[3*j], trianglesB[3*j + 1], trianglesB[3*j + 2]))
                    {
                        found = true;
                        if (a_pairs == NULL)
                        {
                            if (a_stats != NULL)
                            {
                                a_stats->m_boxTests += boxTests;
                                a_stats->m_triangleTests += triangleTests;
                            }
                            return (true);
                        }
                        cOBBTrianglePair hit;
                        hit.m_triangleA = triA;
                        hit.m_triangleB = indicesB[j];
                        a_pairs->push_back(hit);
                    }
                }
            }
            continue;
        }

        // descend into the larger box (or the only internal node)
        const bool descendA = !nodeA.isLeaf() &&
                              (nodeB.isLeaf() || (boxVolume(nodeA.m_box) >= boxVolume(nodeB.m_box)));
        if (descendA)
        {
            stack.push_back(std::make_pair(nodeA.m_right, pair.second));
            stack.push_back(std::make_pair(nodeA.m_left, pair.second));
        }
        else
        {
            stack.push_back(std::make_pair(pair.first, nodeB.m_right));
            stack.push_back(std::make_pair(pair.first, nodeB.m_left));
        }
    }

    if (a_stats != NULL)
    {
        a_stats->m_boxTests += boxTests;
        a_stats->m_triangleTests += triangleTests;
    }

    return (found);
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
