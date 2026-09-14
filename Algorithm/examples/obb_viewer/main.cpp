//==============================================================================
/*
    OBB tree viewer.

    Displays a hand and an object with the wireframes of their OBB trees
    (built with SAH), and colors the hand red when it intersects the object.

    Keys:
        B            show / hide the OBBs
        + / -        next / previous tree level
        L            one level  /  all levels from the root
        C            color by level  /  single color
        T            transparent meshes on / off
        arrows       move the hand in the horizontal plane
        PgUp / PgDn  move the hand up / down
        R            rotate the hand
        Q / Esc      quit
    Mouse:
        left drag    orbit the camera
        wheel        zoom

    Command line:
        --level N            initial tree level
        --cumulative         show all levels from the root
        --hide               start with the OBBs hidden
        --hand x y z         initial position of the hand (meters)
        --triangles N        approximate number of triangles of the object
        --screenshot f.png   save a screenshot and exit
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
#include "CCollisionOBB.h"
//------------------------------------------------------------------------------
#include <GLFW/glfw3.h>
//------------------------------------------------------------------------------
#include <cmath>
#include <iostream>
#include <string>
//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// SCENE
//------------------------------------------------------------------------------

cWorld* world = nullptr;
cCamera* camera = nullptr;
cViewport* viewport = nullptr;
cDirectionalLight* light = nullptr;

cMesh* object = nullptr;
cMesh* hand = nullptr;
cCollisionOBB* objectDetector = nullptr;
cCollisionOBB* handDetector = nullptr;

cFontPtr font;
cLabel* labelStatus = nullptr;
cLabel* labelHelp1 = nullptr;
cLabel* labelHelp2 = nullptr;

cFrequencyCounter freqCounterGraphics;

//------------------------------------------------------------------------------
// WINDOW
//------------------------------------------------------------------------------

GLFWwindow* window = nullptr;
int windowW = 0;
int windowH = 0;
int framebufferW = 0;
int framebufferH = 0;

//------------------------------------------------------------------------------
// STATE
//------------------------------------------------------------------------------

// display of the boxes
bool showBoxes = true;
int displayLevel = 4;
bool cumulative = false;
bool colorByDepth = true;
bool transparent = true;

// pose of the hand (palm facing down, above the ball)
cVector3d handPos(0.0, 0.047, 0.01);
double handYawDeg = 20.0;

// camera orbit
double cameraRadius = 0.30;
double cameraPolarDeg = 60.0;
double cameraAzimuthDeg = 35.0;
bool dragging = false;
double lastCursorX = 0.0;
double lastCursorY = 0.0;

// contact between the hand and the object
bool contact = false;
size_t contactPairs = 0;

// size of the object, construction and query times
int objectTriangles = 18000;
double buildTimeObject = 0.0;
double queryTimeUs = 0.0;

//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

void onWindowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);
void onFrameBufferSizeCallback(GLFWwindow* a_window, int a_width, int a_height);
void onWindowContentScaleCallback(GLFWwindow* a_window, float a_xscale, float a_yscale);
void onErrorCallback(int a_error, const char* a_description);
void onKeyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);
void onMouseButtonCallback(GLFWwindow* a_window, int a_button, int a_action, int a_mods);
void onCursorPosCallback(GLFWwindow* a_window, double a_x, double a_y);
void onScrollCallback(GLFWwindow* a_window, double a_dx, double a_dy);

void renderGraphics(const string& a_screenshot = "");
void applyDisplaySettings();
void updateHandPose();
void updateCamera();


//------------------------------------------------------------------------------
// GEOMETRY
//------------------------------------------------------------------------------

// A simple hand: palm (box) and five fingers (cylinders), in meters.
void createHand(cMesh* a_mesh)
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

// An object of about a_triangles triangles: a tessellated ball resting on a plate.
void createObject(cMesh* a_mesh, int a_triangles)
{
    const unsigned int resolution = (unsigned int)cMax(4.0, floor(sqrt(a_triangles / 2.0) + 0.5));
    cCreateSphere(a_mesh, 0.04, resolution, resolution);
    cCreateBox(a_mesh, 0.14, 0.01, 0.14, cVector3d(0.0, -0.045, 0.0));
}


//==============================================================================

int main(int argc, char* argv[])
{
    //--------------------------------------------------------------------------
    // COMMAND LINE
    //--------------------------------------------------------------------------

    string screenshot;
    for (int i=1; i<argc; i++)
    {
        const string arg = argv[i];
        if ((arg == "--screenshot") && (i + 1 < argc)) screenshot = argv[++i];
        else if ((arg == "--level") && (i + 1 < argc)) displayLevel = atoi(argv[++i]);
        else if (arg == "--cumulative") cumulative = true;
        else if (arg == "--hide") showBoxes = false;
        else if ((arg == "--triangles") && (i + 1 < argc)) objectTriangles = atoi(argv[++i]);
        else if ((arg == "--hand") && (i + 3 < argc))
        {
            handPos.set(atof(argv[i+1]), atof(argv[i+2]), atof(argv[i+3]));
            i += 3;
        }
    }

    cout << endl;
    cout << "-----------------------------------" << endl;
    cout << "OBB tree viewer" << endl;
    cout << "-----------------------------------" << endl;
    cout << "[B] show/hide OBBs  [+/-] level  [L] one/all levels  [C] colors  [T] transparency" << endl;
    cout << "[arrows, PgUp/PgDn] move hand  [R] rotate hand  [mouse] orbit/zoom  [Q] quit" << endl;


    //--------------------------------------------------------------------------
    // OPENGL - WINDOW DISPLAY
    //--------------------------------------------------------------------------

    if (!glfwInit())
    {
        cout << "failed initialization" << endl;
        return (1);
    }
    glfwSetErrorCallback(onErrorCallback);

    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    windowW = (int)(0.8 * mode->height);
    windowH = (int)(0.6 * mode->height);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_SCALE_TO_MONITOR, GLFW_TRUE);

    window = glfwCreateWindow(windowW, windowH, "OBB tree viewer", NULL, NULL);
    if (!window)
    {
        cout << "failed to create window" << endl;
        glfwTerminate();
        return (1);
    }

    glfwSetKeyCallback(window, onKeyCallback);
    glfwSetWindowSizeCallback(window, onWindowSizeCallback);
    glfwSetFramebufferSizeCallback(window, onFrameBufferSizeCallback);
    glfwSetWindowContentScaleCallback(window, onWindowContentScaleCallback);
    glfwSetMouseButtonCallback(window, onMouseButtonCallback);
    glfwSetCursorPosCallback(window, onCursorPosCallback);
    glfwSetScrollCallback(window, onScrollCallback);

    glfwGetWindowSize(window, &windowW, &windowH);
    glfwGetFramebufferSize(window, &framebufferW, &framebufferH);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    if (glewInit() != GLEW_OK)
    {
        cout << "failed to initialize GLEW library" << endl;
        glfwTerminate();
        return (1);
    }


    //--------------------------------------------------------------------------
    // WORLD - CAMERA - LIGHTING
    //--------------------------------------------------------------------------

    world = new cWorld();
    world->m_backgroundColor.setGrayLevel(0.15f);

    camera = new cCamera(world);
    world->addChild(camera);
    camera->setSphericalReferences(cVector3d(0, 0, 0),     // origin
                                   cVector3d(0, 1, 0),     // zenith: y is up
                                   cVector3d(0, 0, 1));    // azimuth reference
    camera->setClippingPlanes(0.01, 10.0);
    camera->setUseMultipassTransparency(true);
    updateCamera();

    light = new cDirectionalLight(world);
    world->addChild(light);
    light->setEnabled(true);
    light->setDir(-0.5, -1.0, -0.7);


    //--------------------------------------------------------------------------
    // OBJECT AND HAND, WITH THEIR OBB TREES
    //--------------------------------------------------------------------------

    object = new cMesh();
    world->addChild(object);
    createObject(object, objectTriangles);
    object->m_material->setBlueCornflower();

    hand = new cMesh();
    world->addChild(hand);
    createHand(hand);

    // build the OBB trees (SAH) and attach them as collision detectors
    cPrecisionClock buildClock;
    buildClock.start(true);
    objectDetector = cCreateOBBCollisionDetector(object);
    buildTimeObject = buildClock.getCurrentTimeSeconds();
    handDetector = cCreateOBBCollisionDetector(hand);
    if ((objectDetector == nullptr) || (handDetector == nullptr))
    {
        cout << "failed to build the OBB trees" << endl;
        return (1);
    }
    objectDetector->m_color.setWhite();
    handDetector->m_color.setYellow();

    cout << "object: " << object->getNumTriangles() << " triangles, "
         << objectDetector->getTree().getStats().m_numNodes << " OBB nodes, depth "
         << objectDetector->getMaxDepth() << endl;
    cout << "hand:   " << hand->getNumTriangles() << " triangles, "
         << handDetector->getTree().getStats().m_numNodes << " OBB nodes, depth "
         << handDetector->getMaxDepth() << endl;

    updateHandPose();
    applyDisplaySettings();


    //--------------------------------------------------------------------------
    // WIDGETS
    //--------------------------------------------------------------------------

    font = NEW_CFONT_CALIBRI_20();

    labelStatus = new cLabel(font);
    camera->m_frontLayer->addChild(labelStatus);

    labelHelp1 = new cLabel(font);
    labelHelp2 = new cLabel(font);
    camera->m_frontLayer->addChild(labelHelp1);
    camera->m_frontLayer->addChild(labelHelp2);
    labelHelp1->setText("[B] show/hide OBBs   [+/-] level   [L] one/all levels   [C] colors   [T] transparency");
    labelHelp2->setText("[arrows, PgUp/PgDn] move hand   [R] rotate hand   [mouse] orbit / zoom   [Q] quit");
    labelHelp1->m_fontColor.setGrayLevel(0.7f);
    labelHelp2->m_fontColor.setGrayLevel(0.7f);


    //--------------------------------------------------------------------------
    // VIEWPORT DISPLAY
    //--------------------------------------------------------------------------

    float contentScaleW, contentScaleH;
    glfwGetWindowContentScale(window, &contentScaleW, &contentScaleH);
    viewport = new cViewport(camera, contentScaleW, contentScaleH);


    //--------------------------------------------------------------------------
    // MAIN GRAPHIC LOOP
    //--------------------------------------------------------------------------

    if (!screenshot.empty())
    {
        // render a few frames so that the window is settled, then save
        for (int i=0; i<10; i++)
        {
            renderGraphics();
            glfwPollEvents();
        }
        renderGraphics(screenshot);
    }
    else
    {
        while (!glfwWindowShouldClose(window))
        {
            renderGraphics();
            glfwPollEvents();
        }
    }

    // note: the viewport is not deleted, as in the CHAI3D examples:
    // cViewport::~cViewport() deletes its view panel a second time (the panel
    // is already deleted with the camera that owns it), which corrupts the heap
    delete world;
    glfwDestroyWindow(window);
    glfwTerminate();

    return (0);
}

//------------------------------------------------------------------------------

void applyDisplaySettings()
{
    // show / hide: standard CHAI3D flag of the meshes
    object->setShowCollisionDetector(showBoxes);
    hand->setShowCollisionDetector(showBoxes);

    // level: positive = one level (cut through the tree), negative = all levels up to it
    objectDetector->setDisplayDepth(cumulative ? -displayLevel : displayLevel);
    handDetector->setDisplayDepth(cumulative ? -displayLevel : displayLevel);

    objectDetector->setColorByDepth(colorByDepth);
    handDetector->setColorByDepth(colorByDepth);

    const float alpha = transparent ? 0.35f : 1.0f;
    object->setUseTransparency(transparent);
    hand->setUseTransparency(transparent);
    object->setTransparencyLevel(alpha);
    hand->setTransparencyLevel(alpha);
}

//------------------------------------------------------------------------------

void updateHandPose()
{
    cMatrix3d palmDown, yaw;
    palmDown.setAxisAngleRotationDeg(cVector3d(1, 0, 0), 90);
    yaw.setAxisAngleRotationDeg(cVector3d(0, 1, 0), handYawDeg);
    hand->setLocalRot(yaw * palmDown);
    hand->setLocalPos(handPos);
}

//------------------------------------------------------------------------------

void updateCamera()
{
    camera->setSphericalDeg(cameraRadius, cameraPolarDeg, cameraAzimuthDeg);
}

//------------------------------------------------------------------------------

void onWindowSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
    windowW = a_width;
    windowH = a_height;
    renderGraphics();
}

//------------------------------------------------------------------------------

void onFrameBufferSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
    framebufferW = a_width;
    framebufferH = a_height;
}

//------------------------------------------------------------------------------

void onWindowContentScaleCallback(GLFWwindow* a_window, float a_xscale, float a_yscale)
{
    if (viewport != nullptr) viewport->setContentScale(a_xscale, a_yscale);
}

//------------------------------------------------------------------------------

void onErrorCallback(int a_error, const char* a_description)
{
    cout << "Error: " << a_description << endl;
}

//------------------------------------------------------------------------------

void onKeyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
    if ((a_action != GLFW_PRESS) && (a_action != GLFW_REPEAT)) return;

    const double step = 0.003;
    const int maxLevel = cMax(objectDetector->getMaxDepth(), handDetector->getMaxDepth());

    switch (a_key)
    {
        case GLFW_KEY_ESCAPE:
        case GLFW_KEY_Q:         glfwSetWindowShouldClose(a_window, GLFW_TRUE); break;

        case GLFW_KEY_B:         showBoxes = !showBoxes; break;
        case GLFW_KEY_EQUAL:
        case GLFW_KEY_KP_ADD:    displayLevel = cMin(displayLevel + 1, maxLevel); break;
        case GLFW_KEY_MINUS:
        case GLFW_KEY_KP_SUBTRACT: displayLevel = cMax(displayLevel - 1, 0); break;
        case GLFW_KEY_L:         cumulative = !cumulative; break;
        case GLFW_KEY_C:         colorByDepth = !colorByDepth; break;
        case GLFW_KEY_T:         transparent = !transparent; break;

        case GLFW_KEY_LEFT:      handPos(0) -= step; break;
        case GLFW_KEY_RIGHT:     handPos(0) += step; break;
        case GLFW_KEY_UP:        handPos(2) -= step; break;
        case GLFW_KEY_DOWN:      handPos(2) += step; break;
        case GLFW_KEY_PAGE_UP:   handPos(1) += step; break;
        case GLFW_KEY_PAGE_DOWN: handPos(1) -= step; break;
        case GLFW_KEY_R:         handYawDeg += 10.0; break;
        default: break;
    }

    applyDisplaySettings();
    updateHandPose();
}

//------------------------------------------------------------------------------

void onMouseButtonCallback(GLFWwindow* a_window, int a_button, int a_action, int a_mods)
{
    if (a_button != GLFW_MOUSE_BUTTON_LEFT) return;
    dragging = (a_action == GLFW_PRESS);
    glfwGetCursorPos(a_window, &lastCursorX, &lastCursorY);
}

//------------------------------------------------------------------------------

void onCursorPosCallback(GLFWwindow* a_window, double a_x, double a_y)
{
    if (!dragging) return;
    cameraAzimuthDeg -= 0.3 * (a_x - lastCursorX);
    cameraPolarDeg = cClamp(cameraPolarDeg - 0.3 * (a_y - lastCursorY), 5.0, 175.0);
    lastCursorX = a_x;
    lastCursorY = a_y;
    updateCamera();
}

//------------------------------------------------------------------------------

void onScrollCallback(GLFWwindow* a_window, double a_dx, double a_dy)
{
    cameraRadius = cClamp(cameraRadius * pow(0.9, a_dy), 0.05, 2.0);
    updateCamera();
}

//------------------------------------------------------------------------------

void renderGraphics(const string& a_screenshot)
{
    if (viewport == nullptr) return;

    /////////////////////////////////////////////////////////////////////
    // HAND - OBJECT INTERSECTION
    /////////////////////////////////////////////////////////////////////

    world->computeGlobalPositions(true);

    std::vector<cOBBTrianglePair> pairs;
    cPrecisionClock queryClock;
    queryClock.start(true);
    contact = cComputeMeshIntersection(hand, object, &pairs);
    queryTimeUs = 1e6 * queryClock.getCurrentTimeSeconds();
    contactPairs = pairs.size();

    if (contact) hand->m_material->setRedCrimson();
    else         hand->m_material->setWhite();
    hand->setTransparencyLevel(transparent ? 0.35f : 1.0f);


    /////////////////////////////////////////////////////////////////////
    // WIDGETS
    /////////////////////////////////////////////////////////////////////

    const int displayH = viewport->getDisplayHeight();

    string status = "object " + cStr((int)object->getNumTriangles()) + " triangles (tree built in " +
                    cStr(buildTimeObject, 2) + " s)";
    status += string("   OBB: ") + (showBoxes ? "ON" : "OFF");
    status += "   level " + cStr(displayLevel) + (cumulative ? " (all levels from root)" : " (one level)");
    status += string("   contact: ") + (contact ? "YES, " + cStr((int)contactPairs) + " triangle pairs" : "no");
    status += " in " + cStr(queryTimeUs, 0) + " us";
    labelStatus->setText(status);
    labelStatus->m_fontColor = contact ? cColorf(1.0f, 0.4f, 0.4f) : cColorf(1.0f, 1.0f, 1.0f);
    labelStatus->setLocalPos(10, displayH - 30);
    labelHelp1->setLocalPos(10, 34);
    labelHelp2->setLocalPos(10, 10);


    /////////////////////////////////////////////////////////////////////
    // RENDER SCENE
    /////////////////////////////////////////////////////////////////////

    world->updateShadowMaps(false, false);
    viewport->renderView(framebufferW, framebufferH);
    glFinish();

    // screenshot: read the whole back buffer (framebuffer pixels, which differ
    // from window coordinates on scaled displays) before it is swapped
    if (!a_screenshot.empty())
    {
        cImagePtr image = cImage::create();
        image->allocate(framebufferW, framebufferH, GL_RGB);
        glPixelStorei(GL_PACK_ALIGNMENT, 1);
        glReadBuffer(GL_BACK);
        glReadPixels(0, 0, framebufferW, framebufferH, GL_RGB, GL_UNSIGNED_BYTE, image->getData());
        if (image->saveToFile(a_screenshot)) cout << "screenshot saved to " << a_screenshot << endl;
        else                                 cout << "failed to save " << a_screenshot << endl;
    }

    GLenum error = glGetError();
    if (error != GL_NO_ERROR) cout << "OpenGL error: " << gluErrorString(error) << endl;

    glfwSwapBuffers(window);
    freqCounterGraphics.signal(1);
}
