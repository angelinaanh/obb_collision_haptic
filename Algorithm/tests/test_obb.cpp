//==============================================================================
/*
    Tests for cCollisionOBBBox.

    Also checks that code in the Algorithm folder can use CHAI3D (cMesh,
    cQuaternion, cMatrix3d, ...) through the "obb" CMake target.
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBBBox.h"
#include "TestUtils.h"
//------------------------------------------------------------------------------
using namespace chai3d;
//------------------------------------------------------------------------------

int main()
{
    // default box: zero size, identity rotation
    cCollisionOBBBox e;
    CHECK(sizeof(cCollisionOBBBox) == 40);
    CHECK(e.m_rotation[3] == 1.0f && e.m_rotation[0] == 0.0f);
    CHECK(approx(e.getAxis(0), cVector3d(1,0,0)));

    // box rotated 90 deg around Z: local X -> world Y
    cQuaternion q; q.fromAxisAngle(cVector3d(0,0,1), C_PI / 2.0);
    cCollisionOBBBox b(cVector3d(1,2,3), cVector3d(2,0.5,1), q);
    CHECK(approx(b.getCenter(), cVector3d(1,2,3)));
    CHECK(approx(b.getHalfExtents(), cVector3d(2,0.5,1)));

    // storage order is (x, y, z, w)
    CHECK(approx(b.m_rotation[2], sin(C_PI/4)) && approx(b.m_rotation[3], cos(C_PI/4)));
    CHECK(approx(b.m_rotation[0], 0) && approx(b.m_rotation[1], 0));
    CHECK(approx(b.getAxis(0), cVector3d(0,1,0)));
    CHECK(approx(b.getAxis(1), cVector3d(-1,0,0)));
    CHECK(approx(b.getAxis(2), cVector3d(0,0,1)));

    // round trip through cQuaternion (w, x, y, z)
    cQuaternion r = b.getRotation();
    CHECK(approx(r.w, q.w) && approx(r.x, q.x) && approx(r.y, q.y) && approx(r.z, q.z));

    // contains: the long axis (half-extent 2) now points along world Y
    CHECK( b.contains(cVector3d(1, 2 + 1.9, 3)));
    CHECK(!b.contains(cVector3d(1 + 1.9, 2, 3)));
    CHECK( b.contains(cVector3d(1 + 0.4, 2, 3 + 0.9)));
    CHECK(!b.contains(cVector3d(1, 2, 3 + 1.1)));

    // construction from a rotation matrix
    cMatrix3d m; q.toRotMat(m);
    cCollisionOBBBox c(cVector3d(0,0,0), cVector3d(1,1,1), m);
    CHECK(approx(c.getAxis(0), cVector3d(0,1,0)));

    // non-normalized input is normalized; degenerate input -> identity
    c.setRotation(0.0f, 0.0f, 2.0f, 2.0f);
    CHECK(approx(c.getRotation().mag(), 1.0));
    c.setRotation(0.0f, 0.0f, 0.0f, 0.0f);
    CHECK(c.m_rotation[3] == 1.0f);

    // box built from the pose of a CHAI3D mesh
    cMesh* mesh = new cMesh();
    mesh->setLocalPos(0.5, -1.0, 2.0);
    mesh->rotateAboutGlobalAxisDeg(cVector3d(0,0,1), 90);
    cCollisionOBBBox fromMesh(mesh->getLocalPos(), cVector3d(1,1,1), mesh->getLocalRot());
    CHECK(approx(fromMesh.getCenter(), cVector3d(0.5, -1.0, 2.0)));
    CHECK(approx(fromMesh.getAxis(0), cVector3d(0,1,0)));
    delete mesh;

    return (testResult());
}
