//==============================================================================
/*
    Tests for the geometry used to draw the OBB wireframes: corners, edges,
    and the nodes selected by the display depth.
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBB.h"
#include "CCollisionOBBDraw.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
#include <set>
#include <utility>
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

int main()
{
    std::mt19937 rng(5);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    std::uniform_real_distribution<double> size(0.01, 2.0);

    //--------------------------------------------------------------------------
    // edge table: 12 distinct edges between corners that differ in one axis
    //--------------------------------------------------------------------------
    {
        std::set<std::pair<int, int> > edges;
        int perAxis[3] = { 0, 0, 0 };
        for (int e=0; e<12; e++)
        {
            const int a = C_OBB_EDGES[e][0];
            const int b = C_OBB_EDGES[e][1];
            const int diff = a ^ b;
            CHECK((diff == 1) || (diff == 2) || (diff == 4));
            if (diff == 1) perAxis[0]++;
            if (diff == 2) perAxis[1]++;
            if (diff == 4) perAxis[2]++;
            edges.insert(std::make_pair(cMin(a, b), cMax(a, b)));
        }
        CHECK(edges.size() == 12);
        CHECK((perAxis[0] == 4) && (perAxis[1] == 4) && (perAxis[2] == 4));
    }

    //--------------------------------------------------------------------------
    // corners of random boxes: position, orientation and size
    //--------------------------------------------------------------------------
    {
        int bad = 0;
        for (int trial=0; trial<1000; trial++)
        {
            const cMatrix3d R = randomRotation(rng);
            const cVector3d half = randomVector(rng, size);
            const cVector3d center = 10.0 * randomVector(rng, uni);
            const cCollisionOBBBox box(center, half, R);

            cVector3d corners[8];
            cGetOBBCorners(box, corners);

            // the stored box (float) is the reference
            const cVector3d c = box.getCenter();
            const cVector3d h = box.getHalfExtents();
            const cMatrix3d rot = box.getRotationMatrix();
            const cVector3d axis[3] = { rot.getCol0(), rot.getCol1(), rot.getCol2() };

            cVector3d mean(0, 0, 0);
            for (int i=0; i<8; i++)
            {
                mean += corners[i] / 8.0;
                for (int k=0; k<3; k++)
                {
                    // corner i is on the + side of axis k if bit k is set
                    const double s = ((i >> k) & 1) ? 1.0 : -1.0;
                    if (!approx(cDot(corners[i] - c, axis[k]), s * h(k), 1e-9)) bad++;
                }
            }
            if (!approx(mean, c, 1e-9)) bad++;

            // every edge has the length of the box along its axis
            for (int e=0; e<12; e++)
            {
                const int a = C_OBB_EDGES[e][0];
                const int b = C_OBB_EDGES[e][1];
                const int k = ((a ^ b) == 1) ? 0 : (((a ^ b) == 2) ? 1 : 2);
                if (!approx(cDistance(corners[a], corners[b]), 2.0 * h(k), 1e-9)) bad++;
            }

            // the stored box approximates the requested one (float storage)
            if (!approx(c, center, 1e-5) || !approx(h, half, 1e-5)) bad++;
        }
        CHECK(bad == 0);
        printf("corners: 1000 random boxes, %d errors\n", bad);
    }

    //--------------------------------------------------------------------------
    // display depth on a real tree
    //--------------------------------------------------------------------------
    {
        cMesh* mesh = new cMesh();
        cCreateSphere(mesh, 0.04, 48, 48);
        cCreateBox(mesh, 0.14, 0.01, 0.14, cVector3d(0.0, -0.045, 0.0));
        cCollisionOBB* detector = cCreateOBBCollisionDetector(mesh);
        CHECK(detector != NULL);

        const cCollisionOBBTree& tree = detector->getTree();
        const std::vector<cCollisionOBBNode>& nodes = tree.getNodes();
        const unsigned int n = tree.getNumTriangles();
        const int maxDepth = detector->getMaxDepth();
        CHECK(maxDepth == tree.getStats().m_maxLeafDepth);

        int badCut = 0, badCumulative = 0;
        for (int d=0; d<=maxDepth + 1; d++)
        {
            // positive depth: the displayed nodes cover every triangle exactly once
            detector->setDisplayDepth(d);
            std::vector<int> cover(n, 0);
            for (size_t i=0; i<nodes.size(); i++)
            {
                if (!detector->isNodeDisplayed(nodes[i])) continue;
                for (unsigned int j=0; j<nodes[i].m_count; j++) cover[tree.getLeafTriangle(nodes[i].m_first + j)]++;
            }
            for (unsigned int t=0; t<n; t++) if (cover[t] != 1) { badCut++; break; }

            // negative depth: all the nodes from the root to depth |d|
            detector->setDisplayDepth(-d);
            int shown = 0, expected = 0;
            for (size_t i=0; i<nodes.size(); i++)
            {
                if (detector->isNodeDisplayed(nodes[i])) shown++;
                if (nodes[i].m_depth <= d) expected++;
            }
            if (shown != expected) badCumulative++;
        }
        CHECK(badCut == 0);
        CHECK(badCumulative == 0);
        printf("display depth: levels 0..%d, %d levels with incomplete cover, %d wrong cumulative levels\n",
               maxDepth + 1, badCut, badCumulative);

        // show / hide goes through the standard CHAI3D flag of the mesh
        mesh->setShowCollisionDetector(true);
        CHECK(mesh->getShowCollisionDetector());
        mesh->setShowCollisionDetector(false);
        CHECK(!mesh->getShowCollisionDetector());

        delete mesh;
    }

    //--------------------------------------------------------------------------
    // level colors
    //--------------------------------------------------------------------------
    {
        const cColorf root = cOBBDepthColor(0.0);
        const cColorf deepest = cOBBDepthColor(1.0);
        CHECK(root.getB() > root.getR());
        CHECK(deepest.getR() > deepest.getB());
        CHECK(approx(cOBBDepthColor(-1.0).getB(), root.getB(), 1e-6));
        CHECK(approx(cOBBDepthColor(2.0).getR(), deepest.getR(), 1e-6));
    }

    return (testResult());
}
