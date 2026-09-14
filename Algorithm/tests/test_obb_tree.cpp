//==============================================================================
/*
    Tests for the OBB tree (SAH construction, segment queries, tree-tree
    queries) and for its integration in CHAI3D (cCollisionOBB).
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBB.h"
#include "CCollisionOBBIntersect.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
#include <algorithm>
#include <cfloat>
#include <string>
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// A simple hand: palm (box) and five fingers (cylinders), in meters.
//------------------------------------------------------------------------------
static void createHand(cMesh* a_mesh, unsigned int a_sides, unsigned int a_segments)
{
    cCreateBox(a_mesh, 0.08, 0.09, 0.02);

    cMatrix3d alongY;
    alongY.setAxisAngleRotationDeg(cVector3d(1, 0, 0), -90);   // cylinder axis z -> y
    const double length[4] = { 0.065, 0.075, 0.07, 0.055 };
    for (int f=0; f<4; f++)
    {
        cCreateCylinder(a_mesh, length[f], 0.008, a_sides, a_segments, 1, true, true,
                        cVector3d(-0.03 + 0.02 * f, 0.045, 0.0), alongY);
    }

    cMatrix3d thumb;
    thumb.setAxisAngleRotationDeg(cVector3d(0, 0, 1), -55);
    cCreateCylinder(a_mesh, 0.05, 0.009, a_sides, a_segments, 1, true, true,
                    cVector3d(0.04, -0.01, 0.0), thumb * alongY);
}

//------------------------------------------------------------------------------
// An object with non-uniform triangle density: a finely tessellated ball
// resting on a coarse plate.
//------------------------------------------------------------------------------
static void createObject(cMesh* a_mesh, unsigned int a_resolution)
{
    cCreateSphere(a_mesh, 0.04, a_resolution, a_resolution);
    cCreateBox(a_mesh, 0.14, 0.01, 0.14, cVector3d(0.0, -0.045, 0.0));
}

//------------------------------------------------------------------------------
// Structural checks: leaves partition the triangles, children split the
// range of their parent, boxes contain their triangles, splits are balanced.
//------------------------------------------------------------------------------
static bool checkTree(const cCollisionOBBTree& a_tree, std::string& a_error)
{
    const std::vector<cCollisionOBBNode>& nodes = a_tree.getNodes();
    const unsigned int n = a_tree.getNumTriangles();
    const cOBBTreeSettings& settings = a_tree.getSettings();

    std::vector<int> seen(n, 0);
    for (size_t i=0; i<nodes.size(); i++)
    {
        const cCollisionOBBNode& node = nodes[i];

        for (unsigned int j=0; j<node.m_count; j++)
        {
            const unsigned int tri = a_tree.getLeafTriangle(node.m_first + j);
            if (node.isLeaf()) seen[tri]++;

            cVector3d v[3];
            a_tree.getTriangle(tri, v[0], v[1], v[2]);
            for (int k=0; k<3; k++)
            {
                if (!node.m_box.contains(v[k])) { a_error = "box does not contain its triangles"; return (false); }
            }
        }

        if (node.isLeaf()) continue;

        const cCollisionOBBNode& L = nodes[node.m_left];
        const cCollisionOBBNode& R = nodes[node.m_right];
        if ((L.m_first != node.m_first) || (R.m_first != node.m_first + L.m_count) ||
            (L.m_count + R.m_count != node.m_count) ||
            (L.m_depth != node.m_depth + 1) || (R.m_depth != node.m_depth + 1))
        {
            a_error = "children do not split the range of their parent";
            return (false);
        }

        if (settings.m_splitMethod == C_OBB_SPLIT_SAH)
        {
            unsigned int minCount = (unsigned int)ceil(settings.m_minSplitFraction * node.m_count);
            minCount = cClamp(minCount, 1u, node.m_count / 2);
            if (cMin(L.m_count, R.m_count) < minCount) { a_error = "unbalanced split"; return (false); }
        }
    }

    for (unsigned int i=0; i<n; i++)
    {
        if (seen[i] != 1) { a_error = "a triangle is not in exactly one leaf"; return (false); }
    }
    return (true);
}

static void printStats(const char* a_name, const cOBBTreeStats& s)
{
    printf("  %-7s nodes %6d, leaves %6d, leaf size mean %.2f max %d, leaf depth %d..%d (mean %.1f),\n"
           "          SAH cost %.1f, node volume sum %.1f x root, sibling overlap %.3f, build %.3f s\n",
           a_name, s.m_numNodes, s.m_numLeaves, s.m_meanLeafSize, s.m_maxLeafSize,
           s.m_minLeafDepth, s.m_maxLeafDepth, s.m_meanLeafDepth,
           s.m_sahCost, s.m_volumeRatio, s.m_meanSiblingOverlap, s.m_buildTime);
}

//------------------------------------------------------------------------------
// Brute force references.
//------------------------------------------------------------------------------
static bool bruteSegment(const cCollisionOBBTree& a_tree, const cVector3d& a, const cVector3d& b, double& a_t)
{
    a_t = DBL_MAX;
    for (unsigned int i=0; i<a_tree.getNumTriangles(); i++)
    {
        cVector3d v0, v1, v2;
        a_tree.getTriangle(i, v0, v1, v2);
        double t;
        if (cTestSegmentTriangle(a, b, v0, v1, v2, t) && (t < a_t)) a_t = t;
    }
    return (a_t < DBL_MAX);
}

static bool pairLess(const cOBBTrianglePair& x, const cOBBTrianglePair& y)
{
    return ((x.m_triangleA < y.m_triangleA) || ((x.m_triangleA == y.m_triangleA) && (x.m_triangleB < y.m_triangleB)));
}

static void brutePairs(const cCollisionOBBTree& A, const cMatrix3d& rotA, const cVector3d& posA,
                       const cCollisionOBBTree& B, const cMatrix3d& rotB, const cVector3d& posB,
                       std::vector<cOBBTrianglePair>& a_pairs)
{
    // same frame convention as cCollisionOBBTree::intersect (B in the frame of A)
    const cMatrix3d rotAT = cTranspose(rotA);
    const cMatrix3d R = rotAT * rotB;
    const cVector3d t = rotAT * (posB - posA);

    std::vector<cVector3d> qB(3 * B.getNumTriangles());
    for (unsigned int j=0; j<B.getNumTriangles(); j++)
    {
        cVector3d q0, q1, q2;
        B.getTriangle(j, q0, q1, q2);
        qB[3*j] = R * q0 + t; qB[3*j+1] = R * q1 + t; qB[3*j+2] = R * q2 + t;
    }

    a_pairs.clear();
    for (unsigned int i=0; i<A.getNumTriangles(); i++)
    {
        cVector3d p0, p1, p2;
        A.getTriangle(i, p0, p1, p2);
        for (unsigned int j=0; j<B.getNumTriangles(); j++)
        {
            if (cTestTriangleTriangle(p0, p1, p2, qB[3*j], qB[3*j+1], qB[3*j+2]))
            {
                cOBBTrianglePair pair = { i, j };
                a_pairs.push_back(pair);
            }
        }
    }
}

static void randomHandPose(std::mt19937& a_rng, cMatrix3d& a_rot, cVector3d& a_pos)
{
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    cVector3d p;
    do { p = randomVector(a_rng, uni); } while (p.length() > 1.0);
    a_pos = 0.09 * p;
    a_rot = randomRotation(a_rng);
}


int main()
{
    std::mt19937 rng(31);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);

    cOBBTreeSettings sah;
    cOBBTreeSettings median;
    median.m_splitMethod = C_OBB_SPLIT_MEDIAN;

    //--------------------------------------------------------------------------
    // construction: structure and statistics
    //--------------------------------------------------------------------------
    cMesh* object = new cMesh();
    createObject(object, 96);
    cMesh* hand = new cMesh();
    createHand(hand, 24, 12);

    cCollisionOBBTree objectSAH, objectMedian, handSAH, handMedian;
    CHECK(objectSAH.build(object, sah));
    CHECK(objectMedian.build(object, median));
    CHECK(handSAH.build(hand, sah));
    CHECK(handMedian.build(hand, median));

    std::string error;
    CHECK(checkTree(objectSAH, error));    if (!error.empty()) printf("  object SAH: %s\n", error.c_str());
    CHECK(checkTree(objectMedian, error)); if (!error.empty()) printf("  object median: %s\n", error.c_str());
    CHECK(checkTree(handSAH, error));      if (!error.empty()) printf("  hand SAH: %s\n", error.c_str());
    CHECK(checkTree(handMedian, error));   if (!error.empty()) printf("  hand median: %s\n", error.c_str());
    CHECK(objectSAH.getStats().m_sahCost < objectMedian.getStats().m_sahCost);

    printf("object (%u triangles):\n", objectSAH.getNumTriangles());
    printStats("SAH", objectSAH.getStats());
    printStats("median", objectMedian.getStats());
    printf("hand (%u triangles):\n", handSAH.getNumTriangles());
    printStats("SAH", handSAH.getStats());
    printStats("median", handMedian.getStats());

    //--------------------------------------------------------------------------
    // segment queries: nearest hit against brute force, cost SAH vs median
    //--------------------------------------------------------------------------
    {
        const int trials = 3000;
        int mismatches = 0, hits = 0;
        cOBBQueryStats statsSAH, statsMedian;
        double timeSAH = 0.0, timeMedian = 0.0;
        cPrecisionClock clock;

        for (int i=0; i<trials; i++)
        {
            const cVector3d a = 0.1 * randomVector(rng, uni);
            const cVector3d b = 0.1 * randomVector(rng, uni);

            cOBBSegmentHit hitSAH, hitMedian;
            clock.start(true);
            objectSAH.intersectSegment(a, b, hitSAH, &statsSAH);
            timeSAH += clock.getCurrentTimeSeconds();
            clock.start(true);
            objectMedian.intersectSegment(a, b, hitMedian, &statsMedian);
            timeMedian += clock.getCurrentTimeSeconds();

            double tBrute;
            const bool brute = bruteSegment(objectSAH, a, b, tBrute);
            if ((hitSAH.m_hit != brute) || (hitMedian.m_hit != brute)) mismatches++;
            else if (brute && (!approx(hitSAH.m_t, tBrute, 1e-12) || !approx(hitMedian.m_t, tBrute, 1e-12))) mismatches++;
            if (brute) hits++;

            // the reported triangle is the one of the mesh
            if (hitSAH.m_hit)
            {
                CHECK(hitSAH.m_meshTriangle == objectSAH.getMeshTriangle(hitSAH.m_triangle));
            }
        }
        CHECK(mismatches == 0);

        const double n = (double)objectSAH.getNumTriangles();
        printf("segment queries (%d, %d hitting), %d mismatches with brute force (%.0f triangle tests each):\n",
               trials, hits, mismatches, n);
        printf("  SAH     %.1f box tests + %.1f triangle tests per query, %.2f us\n",
               (double)statsSAH.m_boxTests / trials, (double)statsSAH.m_triangleTests / trials, 1e6 * timeSAH / trials);
        printf("  median  %.1f box tests + %.1f triangle tests per query, %.2f us\n",
               (double)statsMedian.m_boxTests / trials, (double)statsMedian.m_triangleTests / trials, 1e6 * timeMedian / trials);
    }

    //--------------------------------------------------------------------------
    // CHAI3D integration: same collision events as the AABB detector of CHAI3D
    //--------------------------------------------------------------------------
    {
        cMatrix3d rot;
        rot.setAxisAngleRotationDeg(cVector3d(1, -2, 0.5), 33);
        const cVector3d pos(0.3, -0.1, 0.2);

        const double radii[2] = { 0.0, 0.003 };
        for (int r=0; r<2; r++)
        {
            cMesh* meshOBB = new cMesh();
            cMesh* meshAABB = new cMesh();
            createObject(meshOBB, 48);
            createObject(meshAABB, 48);
            meshOBB->setLocalPos(pos);   meshOBB->setLocalRot(rot);
            meshAABB->setLocalPos(pos);  meshAABB->setLocalRot(rot);

            CHECK(cCreateOBBCollisionDetector(meshOBB) != NULL);
            CHECK(dynamic_cast<cCollisionOBB*>(meshOBB->getCollisionDetector()) != NULL);
            meshAABB->createAABBCollisionDetector(radii[r]);

            cCollisionSettings settings;
            settings.m_collisionRadius = radii[r];
            settings.m_checkForNearestCollisionOnly = false;

            int mismatches = 0, hits = 0;
            const int trials = 2000;
            for (int i=0; i<trials; i++)
            {
                const cVector3d a = pos + 0.1 * randomVector(rng, uni);
                const cVector3d b = pos + 0.1 * randomVector(rng, uni);

                cCollisionRecorder recOBB, recAABB;
                const bool hitOBB = meshOBB->computeCollisionDetection(a, b, recOBB, settings);
                const bool hitAABB = meshAABB->computeCollisionDetection(a, b, recAABB, settings);

                std::vector<int> idOBB, idAABB;
                for (size_t k=0; k<recOBB.m_collisions.size(); k++) idOBB.push_back(recOBB.m_collisions[k].m_index);
                for (size_t k=0; k<recAABB.m_collisions.size(); k++) idAABB.push_back(recAABB.m_collisions[k].m_index);
                std::sort(idOBB.begin(), idOBB.end());
                std::sort(idAABB.begin(), idAABB.end());

                if ((hitOBB != hitAABB) || (idOBB != idAABB) ||
                    (hitOBB && !approx(recOBB.m_nearestCollision.m_squareDistance,
                                       recAABB.m_nearestCollision.m_squareDistance, 1e-15)))
                {
                    mismatches++;
                }
                if (hitAABB) hits++;
            }
            CHECK(mismatches == 0);
            printf("CHAI3D detector, tool radius %.3f: %d segments (%d hitting), %d differences with cCollisionAABB\n",
                   radii[r], trials, hits, mismatches);

            delete meshOBB;
            delete meshAABB;
        }
    }

    //--------------------------------------------------------------------------
    // hand vs object: all intersecting triangle pairs against brute force
    //--------------------------------------------------------------------------
    {
        cMesh* smallObject = new cMesh();
        createObject(smallObject, 24);
        cMesh* smallHand = new cMesh();
        createHand(smallHand, 12, 4);

        cCollisionOBBTree treeObject, treeHand;
        treeObject.build(smallObject, sah);
        treeHand.build(smallHand, sah);

        const int poses = 25;
        int mismatches = 0, touching = 0;
        cOBBQueryStats stats;
        for (int i=0; i<poses; i++)
        {
            cMatrix3d rot;
            cVector3d pos;
            randomHandPose(rng, rot, pos);

            std::vector<cOBBTrianglePair> pairs, reference;
            const bool hit = cCollisionOBBTree::intersect(treeHand, rot, pos, treeObject, cIdentity3d(), cVector3d(0, 0, 0), &pairs, &stats);
            const bool first = cCollisionOBBTree::intersect(treeHand, rot, pos, treeObject, cIdentity3d(), cVector3d(0, 0, 0));
            brutePairs(treeHand, rot, pos, treeObject, cIdentity3d(), cVector3d(0, 0, 0), reference);

            std::sort(pairs.begin(), pairs.end(), pairLess);
            std::sort(reference.begin(), reference.end(), pairLess);
            bool same = (pairs.size() == reference.size());
            for (size_t k=0; same && k<pairs.size(); k++)
            {
                same = (pairs[k].m_triangleA == reference[k].m_triangleA) && (pairs[k].m_triangleB == reference[k].m_triangleB);
            }
            if (!same || (hit != !reference.empty()) || (first != hit)) mismatches++;
            if (hit) touching++;
        }
        CHECK(mismatches == 0);
        CHECK(touching > 0 && touching < poses);
        printf("hand (%u tri) vs object (%u tri): %d poses (%d in contact), %d mismatches with brute force,\n"
               "  %.0f box tests + %.0f triangle tests per query instead of %u triangle tests\n",
               treeHand.getNumTriangles(), treeObject.getNumTriangles(), poses, touching, mismatches,
               (double)stats.m_boxTests / poses, (double)stats.m_triangleTests / poses,
               treeHand.getNumTriangles() * treeObject.getNumTriangles());

        delete smallObject;
        delete smallHand;
    }

    //--------------------------------------------------------------------------
    // hand vs object through the meshes (global poses), and query cost
    //--------------------------------------------------------------------------
    {
        CHECK(cCreateOBBCollisionDetector(object) != NULL);
        CHECK(cCreateOBBCollisionDetector(hand) != NULL);

        const int poses = 500;
        cOBBQueryStats statsFirst, statsAll, statsMedian;
        double timeFirst = 0.0, timeAll = 0.0, timeMedian = 0.0;
        int contacts = 0, consistent = 0;
        long long pairsTotal = 0;
        cPrecisionClock clock;

        for (int i=0; i<poses; i++)
        {
            cMatrix3d rot;
            cVector3d pos;
            randomHandPose(rng, rot, pos);
            hand->setLocalRot(rot);
            hand->setLocalPos(pos);
            hand->computeGlobalPositions(true);
            object->computeGlobalPositions(true);

            std::vector<cOBBTrianglePair> pairs, pairsMedian;
            clock.start(true);
            const bool first = cComputeMeshIntersection(hand, object, NULL, &statsFirst);
            timeFirst += clock.getCurrentTimeSeconds();
            clock.start(true);
            const bool all = cComputeMeshIntersection(hand, object, &pairs, &statsAll);
            timeAll += clock.getCurrentTimeSeconds();
            clock.start(true);
            cCollisionOBBTree::intersect(handMedian, rot, pos, objectMedian, cIdentity3d(), cVector3d(0, 0, 0), &pairsMedian, &statsMedian);
            timeMedian += clock.getCurrentTimeSeconds();

            // pairs are reported with mesh triangle indices
            bool valid = true;
            for (size_t k=0; k<pairs.size(); k++)
            {
                valid = valid && (pairs[k].m_triangleA < hand->getNumTriangles()) && (pairs[k].m_triangleB < object->getNumTriangles());
            }
            if ((first == all) && (pairs.size() == pairsMedian.size()) && valid) consistent++;
            if (all) contacts++;
            pairsTotal += (long long)pairs.size();
        }
        CHECK(consistent == poses);
        printf("hand (%u tri) vs object (%u tri) via cComputeMeshIntersection, %d poses (%d in contact, %.1f pairs on average):\n",
               hand->getNumTriangles(), object->getNumTriangles(), poses, contacts, (double)pairsTotal / cMax(contacts, 1));
        printf("  first contact, SAH  %8.1f box tests + %7.1f triangle tests, %8.1f us\n",
               (double)statsFirst.m_boxTests / poses, (double)statsFirst.m_triangleTests / poses, 1e6 * timeFirst / poses);
        printf("  all pairs,     SAH  %8.1f box tests + %7.1f triangle tests, %8.1f us\n",
               (double)statsAll.m_boxTests / poses, (double)statsAll.m_triangleTests / poses, 1e6 * timeAll / poses);
        printf("  all pairs,  median  %8.1f box tests + %7.1f triangle tests, %8.1f us\n",
               (double)statsMedian.m_boxTests / poses, (double)statsMedian.m_triangleTests / poses, 1e6 * timeMedian / poses);
    }

    delete object;
    delete hand;

    return (testResult());
}
