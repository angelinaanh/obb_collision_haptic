//==============================================================================
/*
    OBB tree benchmark.

    For objects of about 1,000 / 10,000 / 100,000 / 1,000,000 triangles, measures:
      - construction time: OBB tree (SAH), OBB tree (median split), CHAI3D AABB;
      - segment queries (haptic tool): OBB tree vs CHAI3D AABB vs brute force,
        and agreement with brute force;
      - hand vs object queries: first contact, all pairs, and agreement with an
        exact reference that uses no OBB and no tree (uniform grid + exact
        triangle-triangle tests).

    Usage: obb_benchmark [results.csv] [maxTriangles]
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBB.h"
#include "CCollisionOBBIntersect.h"
//------------------------------------------------------------------------------
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>
//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// SCENE
//------------------------------------------------------------------------------

// A simple hand: palm (box) and five fingers (cylinders), in meters.
static void createHand(cMesh* a_mesh)
{
    cCreateBox(a_mesh, 0.08, 0.09, 0.02);

    cMatrix3d alongY;
    alongY.setAxisAngleRotationDeg(cVector3d(1, 0, 0), -90);
    const double length[4] = { 0.065, 0.075, 0.07, 0.055 };
    for (int f=0; f<4; f++)
    {
        cCreateCylinder(a_mesh, length[f], 0.008, 24, 12, 1, true, true,
                        cVector3d(-0.03 + 0.02 * f, 0.045, 0.0), alongY);
    }

    cMatrix3d thumb;
    thumb.setAxisAngleRotationDeg(cVector3d(0, 0, 1), -55);
    cCreateCylinder(a_mesh, 0.05, 0.009, 24, 12, 1, true, true,
                    cVector3d(0.04, -0.01, 0.0), thumb * alongY);
}

// An object of about a_triangles triangles: a tessellated ball on a plate.
static void createObject(cMesh* a_mesh, int a_triangles)
{
    const unsigned int resolution = (unsigned int)cMax(4.0, floor(sqrt(a_triangles / 2.0) + 0.5));
    cCreateSphere(a_mesh, 0.04, resolution, resolution);
    cCreateBox(a_mesh, 0.14, 0.01, 0.14, cVector3d(0.0, -0.045, 0.0));
    a_mesh->m_material->setHapticTriangleSides(true, true);
}

//------------------------------------------------------------------------------
// RANDOM (draws sequenced explicitly, see tests/TestUtils.h)
//------------------------------------------------------------------------------

static cVector3d randomVector(mt19937& a_rng, uniform_real_distribution<double>& a_dist)
{
    const double x = a_dist(a_rng);
    const double y = a_dist(a_rng);
    const double z = a_dist(a_rng);
    return (cVector3d(x, y, z));
}

static cMatrix3d randomRotation(mt19937& a_rng)
{
    normal_distribution<double> gauss(0.0, 1.0);
    const double w = gauss(a_rng), x = gauss(a_rng), y = gauss(a_rng), z = gauss(a_rng);
    cQuaternion q(w, x, y, z);
    q.normalize();
    cMatrix3d rot;
    q.toRotMat(rot);
    return (rot);
}

static void randomHandPose(mt19937& a_rng, cMatrix3d& a_rot, cVector3d& a_pos)
{
    uniform_real_distribution<double> uni(-1.0, 1.0);
    cVector3d p;
    do { p = randomVector(a_rng, uni); } while (p.length() > 1.0);
    a_pos = 0.09 * p;
    a_rot = randomRotation(a_rng);
}

//------------------------------------------------------------------------------
// EXACT REFERENCE FOR HAND vs OBJECT (no OBB, no tree)
//------------------------------------------------------------------------------

struct cBounds
{
    cVector3d m_lo, m_hi;
    cBounds() : m_lo(DBL_MAX, DBL_MAX, DBL_MAX), m_hi(-DBL_MAX, -DBL_MAX, -DBL_MAX) {}
    void enclose(const cVector3d& p)
    {
        for (int k=0; k<3; k++) { m_lo(k) = cMin(m_lo(k), p(k)); m_hi(k) = cMax(m_hi(k), p(k)); }
    }
    bool overlaps(const cBounds& b) const
    {
        for (int k=0; k<3; k++) if ((m_hi(k) < b.m_lo(k)) || (b.m_hi(k) < m_lo(k))) return (false);
        return (true);
    }
};

// Uniform grid over the hand triangles, in the local frame of the hand.
struct cHandGrid
{
    static const int N = 24;
    cBounds m_bounds;
    double m_cell[3];
    vector<cVector3d> m_tri;        // 3 vertices per hand triangle
    vector<cBounds> m_triBounds;
    vector<vector<unsigned int> > m_cells;

    void cellRange(const cBounds& b, int lo[3], int hi[3]) const
    {
        for (int k=0; k<3; k++)
        {
            lo[k] = cClamp((int)floor((b.m_lo(k) - m_bounds.m_lo(k)) / m_cell[k]), 0, N - 1);
            hi[k] = cClamp((int)floor((b.m_hi(k) - m_bounds.m_lo(k)) / m_cell[k]), 0, N - 1);
        }
    }

    void build(const cCollisionOBBTree& a_hand)
    {
        const unsigned int n = a_hand.getNumTriangles();
        m_tri.resize(3 * n);
        m_triBounds.resize(n);
        for (unsigned int i=0; i<n; i++)
        {
            a_hand.getTriangle(i, m_tri[3*i], m_tri[3*i+1], m_tri[3*i+2]);
            for (int k=0; k<3; k++) { m_triBounds[i].enclose(m_tri[3*i+k]); m_bounds.enclose(m_tri[3*i+k]); }
        }
        for (int k=0; k<3; k++) m_cell[k] = cMax((m_bounds.m_hi(k) - m_bounds.m_lo(k)) / N, 1e-12);

        m_cells.assign(N * N * N, vector<unsigned int>());
        for (unsigned int i=0; i<n; i++)
        {
            int lo[3], hi[3];
            cellRange(m_triBounds[i], lo, hi);
            for (int x=lo[0]; x<=hi[0]; x++)
                for (int y=lo[1]; y<=hi[1]; y++)
                    for (int z=lo[2]; z<=hi[2]; z++)
                        m_cells[(x * N + y) * N + z].push_back(i);
        }
    }
};

// All intersecting pairs (hand triangle, object triangle) for the object placed
// in the frame of the hand with the same arithmetic as cCollisionOBBTree::intersect().
static void referencePairs(const cHandGrid& a_grid,
                           const cCollisionOBBTree& a_object,
                           const vector<cBounds>& a_objectBoundsWorld,
                           const cMatrix3d& a_rotHand, const cVector3d& a_posHand,
                           vector<cOBBTrianglePair>& a_pairs)
{
    const cMatrix3d rotAT = cTranspose(a_rotHand);
    const cMatrix3d R = rotAT * cIdentity3d();
    const cVector3d t = rotAT * (cVector3d(0, 0, 0) - a_posHand);

    // bounds of the hand in world coordinates (conservative)
    cBounds handWorld;
    for (int c=0; c<8; c++)
    {
        const cVector3d corner((c & 1) ? a_grid.m_bounds.m_hi(0) : a_grid.m_bounds.m_lo(0),
                               (c & 2) ? a_grid.m_bounds.m_hi(1) : a_grid.m_bounds.m_lo(1),
                               (c & 4) ? a_grid.m_bounds.m_hi(2) : a_grid.m_bounds.m_lo(2));
        handWorld.enclose(a_rotHand * corner + a_posHand);
    }

    a_pairs.clear();
    vector<unsigned int> stamp(a_grid.m_triBounds.size(), 0);
    unsigned int current = 0;

    for (unsigned int j=0; j<a_object.getNumTriangles(); j++)
    {
        if (!a_objectBoundsWorld[j].overlaps(handWorld)) continue;

        cVector3d q[3];
        a_object.getTriangle(j, q[0], q[1], q[2]);
        cBounds qb;
        for (int k=0; k<3; k++) { q[k] = R * q[k] + t; qb.enclose(q[k]); }
        if (!qb.overlaps(a_grid.m_bounds)) continue;

        current++;
        int lo[3], hi[3];
        a_grid.cellRange(qb, lo, hi);
        for (int x=lo[0]; x<=hi[0]; x++)
            for (int y=lo[1]; y<=hi[1]; y++)
                for (int z=lo[2]; z<=hi[2]; z++)
                {
                    const vector<unsigned int>& cell = a_grid.m_cells[(x * cHandGrid::N + y) * cHandGrid::N + z];
                    for (size_t c=0; c<cell.size(); c++)
                    {
                        const unsigned int i = cell[c];
                        if (stamp[i] == current) continue;
                        stamp[i] = current;
                        if (!a_grid.m_triBounds[i].overlaps(qb)) continue;
                        if (cTestTriangleTriangle(a_grid.m_tri[3*i], a_grid.m_tri[3*i+1], a_grid.m_tri[3*i+2], q[0], q[1], q[2]))
                        {
                            cOBBTrianglePair pair = { i, j };
                            a_pairs.push_back(pair);
                        }
                    }
                }
    }
}

static bool pairLess(const cOBBTrianglePair& x, const cOBBTrianglePair& y)
{
    return ((x.m_triangleA < y.m_triangleA) || ((x.m_triangleA == y.m_triangleA) && (x.m_triangleB < y.m_triangleB)));
}

//------------------------------------------------------------------------------
// BENCHMARK
//------------------------------------------------------------------------------

struct cRow
{
    int n = 0, handTriangles = 0;
    double buildSAH = 0, buildMedian = 0, buildAABB = 0;
    int nodesSAH = 0, depthSAH = 0, depthMedian = 0;
    double sahCostSAH = 0, sahCostMedian = 0, memSAH = 0, memAABB = 0;

    int segments = 0, segHits = 0, segRef = 0;
    double segOBBChai3d = 0, segAABBChai3d = 0, segSAH = 0, segMedian = 0, segBrute = 0;
    double segBoxSAH = 0, segTriSAH = 0, segBoxMedian = 0, segTriMedian = 0;
    int segWrongSAH = 0, segWrongOBBChai3d = 0, segWrongAABBChai3d = 0, segDiffOBBvsAABB = 0;

    // detection on the segments checked against brute force: true hits found
    // (same nearest point) and hits reported, per method
    int segBruteHits = 0;
    int segFound[3] = { 0, 0, 0 };      // OBB direct, OBB via CHAI3D, CHAI3D AABB
    int segReported[3] = { 0, 0, 0 };

    int poses = 0, contacts = 0;
    double meanPairs = 0, handFirst = 0, handAll = 0, handAllMedian = 0, handRef = 0;
    double handBoxAll = 0, handTriAll = 0, precision = 0, recall = 0;
    int contactDisagree = 0;
};

static cRow runSize(int a_target, mt19937& a_rng, cMesh* a_hand, const cCollisionOBBTree& a_handSAH,
                    const cCollisionOBBTree& a_handMedian, const cHandGrid& a_grid)
{
    cRow row;
    cPrecisionClock clock;
    uniform_real_distribution<double> uni(-1.0, 1.0);

    cOBBTreeSettings median;
    median.m_splitMethod = C_OBB_SPLIT_MEDIAN;

    //--------------------------------------------------------------------------
    // construction
    //--------------------------------------------------------------------------
    cMesh* meshOBB = new cMesh();
    cMesh* meshAABB = new cMesh();
    createObject(meshOBB, a_target);
    createObject(meshAABB, a_target);
    row.n = (int)meshOBB->getNumTriangles();
    row.handTriangles = (int)a_hand->getNumTriangles();
    printf("\n=== object: %d triangles ===\n", row.n);

    clock.start(true);
    cCollisionOBB* detector = cCreateOBBCollisionDetector(meshOBB);
    row.buildSAH = clock.getCurrentTimeSeconds();
    const cCollisionOBBTree& tree = detector->getTree();

    clock.start(true);
    meshAABB->createAABBCollisionDetector(0.0);
    row.buildAABB = clock.getCurrentTimeSeconds();

    cCollisionOBBTree treeMedian;
    clock.start(true);
    treeMedian.build(meshOBB, median);
    row.buildMedian = clock.getCurrentTimeSeconds();

    const cOBBTreeStats& s = tree.getStats();
    row.nodesSAH = s.m_numNodes;
    row.depthSAH = s.m_maxLeafDepth;
    row.depthMedian = treeMedian.getStats().m_maxLeafDepth;
    row.sahCostSAH = s.m_sahCost;
    row.sahCostMedian = treeMedian.getStats().m_sahCost;
    row.memSAH = (s.m_numNodes * sizeof(cCollisionOBBNode) + row.n * (3 * sizeof(unsigned int) + 2 * sizeof(int)) +
                  meshOBB->getNumVertices() * sizeof(cVector3d)) / 1048576.0;
    row.memAABB = ((2.0 * row.n - 1.0) * sizeof(cCollisionAABBNode)) / 1048576.0;
    printf("build: SAH %.3f s (%d nodes, depth %d), median %.3f s, CHAI3D AABB %.3f s\n",
           row.buildSAH, row.nodesSAH, row.depthSAH, row.buildMedian, row.buildAABB);

    //--------------------------------------------------------------------------
    // segment queries
    //--------------------------------------------------------------------------
    const int K = 1000;
    const int refBudget = (int)cClamp(2.0e8 / row.n, 50.0, (double)K);
    row.segments = K;
    row.segRef = refBudget;

    cCollisionSettings settings;
    cOBBQueryStats statsSAH, statsMedian;
    for (int i=0; i<K; i++)
    {
        const cVector3d a = 0.1 * randomVector(a_rng, uni);
        const cVector3d b = 0.1 * randomVector(a_rng, uni);

        cCollisionRecorder recOBB, recAABB;
        clock.start(true);
        const bool hitOBB = meshOBB->computeCollisionDetection(a, b, recOBB, settings);
        row.segOBBChai3d += clock.getCurrentTimeSeconds();
        clock.start(true);
        const bool hitAABB = meshAABB->computeCollisionDetection(a, b, recAABB, settings);
        row.segAABBChai3d += clock.getCurrentTimeSeconds();

        cOBBSegmentHit hitSAH, hitMedian;
        clock.start(true);
        tree.intersectSegment(a, b, hitSAH, &statsSAH);
        row.segSAH += clock.getCurrentTimeSeconds();
        clock.start(true);
        treeMedian.intersectSegment(a, b, hitMedian, &statsMedian);
        row.segMedian += clock.getCurrentTimeSeconds();

        if ((hitOBB != hitAABB) ||
            (hitOBB && fabs(recOBB.m_nearestCollision.m_squareDistance - recAABB.m_nearestCollision.m_squareDistance) > 1e-15))
        {
            row.segDiffOBBvsAABB++;
        }
        if (hitSAH.m_hit) row.segHits++;

        // brute force reference on the first segments
        if (i < refBudget)
        {
            clock.start(true);
            double tBrute = DBL_MAX;
            for (unsigned int j=0; j<tree.getNumTriangles(); j++)
            {
                cVector3d v0, v1, v2;
                tree.getTriangle(j, v0, v1, v2);
                double t;
                if (cTestSegmentTriangle(a, b, v0, v1, v2, t) && (t < tBrute)) tBrute = t;
            }
            row.segBrute += clock.getCurrentTimeSeconds();

            const bool brute = (tBrute < DBL_MAX);
            const double d2 = brute ? cSqr(tBrute * cDistance(a, b)) : 0.0;
            const bool okSAH = (hitSAH.m_hit == brute) && (!brute || fabs(hitSAH.m_t - tBrute) <= 1e-12);
            const bool okOBB = (hitOBB == brute) && (!brute || fabs(recOBB.m_nearestCollision.m_squareDistance - d2) <= 1e-9 * cMax(d2, 1e-12));
            const bool okAABB = (hitAABB == brute) && (!brute || fabs(recAABB.m_nearestCollision.m_squareDistance - d2) <= 1e-9 * cMax(d2, 1e-12));
            if (!okSAH) row.segWrongSAH++;
            if (!okOBB) row.segWrongOBBChai3d++;
            if (!okAABB) row.segWrongAABBChai3d++;

            if (brute) row.segBruteHits++;
            if (hitSAH.m_hit) row.segReported[0]++;
            if (hitOBB) row.segReported[1]++;
            if (hitAABB) row.segReported[2]++;
            if (brute && okSAH) row.segFound[0]++;
            if (brute && okOBB) row.segFound[1]++;
            if (brute && okAABB) row.segFound[2]++;
        }
    }
    row.segOBBChai3d *= 1e6 / K;
    row.segAABBChai3d *= 1e6 / K;
    row.segSAH *= 1e6 / K;
    row.segMedian *= 1e6 / K;
    row.segBrute *= 1e6 / refBudget;
    row.segBoxSAH = (double)statsSAH.m_boxTests / K;
    row.segTriSAH = (double)statsSAH.m_triangleTests / K;
    row.segBoxMedian = (double)statsMedian.m_boxTests / K;
    row.segTriMedian = (double)statsMedian.m_triangleTests / K;
    printf("segments: OBB (CHAI3D API) %.2f us, AABB (CHAI3D) %.2f us, OBB direct %.2f us, brute force %.1f us\n",
           row.segOBBChai3d, row.segAABBChai3d, row.segSAH, row.segBrute);
    printf("          errors vs brute force (%d segments): OBB %d, OBB via CHAI3D %d, AABB %d; OBB vs AABB differences %d/%d\n",
           refBudget, row.segWrongSAH, row.segWrongOBBChai3d, row.segWrongAABBChai3d, row.segDiffOBBvsAABB, K);
    printf("          hits found / true hits (%d): OBB %d, OBB via CHAI3D %d, AABB %d; reported: %d, %d, %d\n",
           row.segBruteHits, row.segFound[0], row.segFound[1], row.segFound[2],
           row.segReported[0], row.segReported[1], row.segReported[2]);

    //--------------------------------------------------------------------------
    // hand vs object
    //--------------------------------------------------------------------------
    vector<cBounds> objectBounds(tree.getNumTriangles());
    for (unsigned int j=0; j<tree.getNumTriangles(); j++)
    {
        cVector3d v0, v1, v2;
        tree.getTriangle(j, v0, v1, v2);
        objectBounds[j].enclose(v0); objectBounds[j].enclose(v1); objectBounds[j].enclose(v2);
    }

    const int P = 200;
    row.poses = P;
    long long found = 0, expected = 0, common = 0, totalPairs = 0;
    cOBBQueryStats statsAll;
    for (int i=0; i<P; i++)
    {
        cMatrix3d rot;
        cVector3d pos;
        randomHandPose(a_rng, rot, pos);

        clock.start(true);
        const bool first = cCollisionOBBTree::intersect(a_handSAH, rot, pos, tree, cIdentity3d(), cVector3d(0, 0, 0));
        row.handFirst += clock.getCurrentTimeSeconds();

        vector<cOBBTrianglePair> pairs, pairsMedian, reference;
        clock.start(true);
        cCollisionOBBTree::intersect(a_handSAH, rot, pos, tree, cIdentity3d(), cVector3d(0, 0, 0), &pairs, &statsAll);
        row.handAll += clock.getCurrentTimeSeconds();

        clock.start(true);
        cCollisionOBBTree::intersect(a_handMedian, rot, pos, treeMedian, cIdentity3d(), cVector3d(0, 0, 0), &pairsMedian);
        row.handAllMedian += clock.getCurrentTimeSeconds();

        clock.start(true);
        referencePairs(a_grid, tree, objectBounds, rot, pos, reference);
        row.handRef += clock.getCurrentTimeSeconds();

        sort(pairs.begin(), pairs.end(), pairLess);
        sort(reference.begin(), reference.end(), pairLess);
        size_t p = 0, q = 0;
        while ((p < pairs.size()) && (q < reference.size()))
        {
            if (pairLess(pairs[p], reference[q])) p++;
            else if (pairLess(reference[q], pairs[p])) q++;
            else { common++; p++; q++; }
        }
        found += (long long)pairs.size();
        expected += (long long)reference.size();
        if (first != !reference.empty()) row.contactDisagree++;
        if (!reference.empty()) { row.contacts++; totalPairs += (long long)reference.size(); }
    }
    row.handFirst *= 1e6 / P;
    row.handAll *= 1e6 / P;
    row.handAllMedian *= 1e6 / P;
    row.handRef *= 1e6 / P;
    row.handBoxAll = (double)statsAll.m_boxTests / P;
    row.handTriAll = (double)statsAll.m_triangleTests / P;
    row.meanPairs = (row.contacts > 0) ? (double)totalPairs / row.contacts : 0.0;
    row.precision = (found > 0) ? (double)common / found : 1.0;
    row.recall = (expected > 0) ? (double)common / expected : 1.0;
    printf("hand vs object (%d poses, %d in contact, %.0f pairs): first contact %.1f us, all pairs %.1f us (median tree %.1f us), reference %.1f us\n",
           P, row.contacts, row.meanPairs, row.handFirst, row.handAll, row.handAllMedian, row.handRef);
    printf("          precision %.6f, recall %.6f, contact disagreements %d\n", row.precision, row.recall, row.contactDisagree);

    delete meshOBB;
    delete meshAABB;
    return (row);
}


int main(int argc, char* argv[])
{
    const string csv = (argc > 1) ? argv[1] : "obb_benchmark.csv";
    const int maxTriangles = (argc > 2) ? atoi(argv[2]) : 1000000;

    mt19937 rng(2026);

    cMesh* hand = new cMesh();
    createHand(hand);
    cCollisionOBBTree handSAH, handMedian;
    cOBBTreeSettings median;
    median.m_splitMethod = C_OBB_SPLIT_MEDIAN;
    handSAH.build(hand);
    handMedian.build(hand, median);
    cHandGrid grid;
    grid.build(handSAH);
    printf("hand: %u triangles\n", hand->getNumTriangles());

    vector<cRow> rows;
    const int sizes[4] = { 1000, 10000, 100000, 1000000 };
    for (int i=0; i<4; i++)
    {
        if (sizes[i] > maxTriangles) break;
        rows.push_back(runSize(sizes[i], rng, hand, handSAH, handMedian, grid));
    }

    ofstream out(csv.c_str());
    out << "triangles,hand_triangles,build_sah_s,build_median_s,build_aabb_s,nodes_sah,depth_sah,depth_median,"
           "sah_cost_sah,sah_cost_median,mem_obb_mb,mem_aabb_mb,"
           "segments,segment_hits,segments_checked,seg_obb_chai3d_us,seg_aabb_chai3d_us,seg_obb_direct_us,seg_median_direct_us,seg_brute_us,"
           "seg_box_tests_sah,seg_tri_tests_sah,seg_box_tests_median,seg_tri_tests_median,"
           "seg_errors_obb,seg_errors_obb_chai3d,seg_errors_aabb,seg_diff_obb_aabb,"
           "poses,contacts,mean_pairs,hand_first_us,hand_all_us,hand_all_median_us,hand_reference_us,"
           "hand_box_tests,hand_tri_tests,precision,recall,contact_disagreements,"
           "seg_true_hits,seg_found_obb,seg_found_obb_chai3d,seg_found_aabb,seg_reported_obb,seg_reported_obb_chai3d,seg_reported_aabb\n";
    for (size_t i=0; i<rows.size(); i++)
    {
        const cRow& r = rows[i];
        char line[2048];
        snprintf(line, sizeof(line),
                 "%d,%d,%.6f,%.6f,%.6f,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,"
                 "%d,%d,%d,%.4f,%.4f,%.4f,%.4f,%.3f,%.2f,%.2f,%.2f,%.2f,%d,%d,%d,%d,"
                 "%d,%d,%.2f,%.3f,%.3f,%.3f,%.3f,%.1f,%.1f,%.6f,%.6f,%d,"
                 "%d,%d,%d,%d,%d,%d,%d\n",
                 r.n, r.handTriangles, r.buildSAH, r.buildMedian, r.buildAABB, r.nodesSAH, r.depthSAH, r.depthMedian,
                 r.sahCostSAH, r.sahCostMedian, r.memSAH, r.memAABB,
                 r.segments, r.segHits, r.segRef, r.segOBBChai3d, r.segAABBChai3d, r.segSAH, r.segMedian, r.segBrute,
                 r.segBoxSAH, r.segTriSAH, r.segBoxMedian, r.segTriMedian,
                 r.segWrongSAH, r.segWrongOBBChai3d, r.segWrongAABBChai3d, r.segDiffOBBvsAABB,
                 r.poses, r.contacts, r.meanPairs, r.handFirst, r.handAll, r.handAllMedian, r.handRef,
                 r.handBoxAll, r.handTriAll, r.precision, r.recall, r.contactDisagree,
                 r.segBruteHits, r.segFound[0], r.segFound[1], r.segFound[2],
                 r.segReported[0], r.segReported[1], r.segReported[2]);
        out << line;
    }
    printf("\nresults written to %s\n", csv.c_str());

    delete hand;
    return (0);
}
