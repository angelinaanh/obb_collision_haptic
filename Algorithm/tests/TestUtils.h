//==============================================================================
/*
    Minimal helpers shared by the OBB tests.
*/
//==============================================================================

//------------------------------------------------------------------------------
#ifndef TestUtilsH
#define TestUtilsH
//------------------------------------------------------------------------------
#include "math/CVector3d.h"
#include "math/CMatrix3d.h"
#include "math/CQuaternion.h"
//------------------------------------------------------------------------------
#include <cmath>
#include <cstdio>
#include <random>
//------------------------------------------------------------------------------

static int g_fails = 0;

#define CHECK(c) do { if (!(c)) { printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); g_fails++; } } while (0)

// note: "near" cannot be used as a name on Windows (empty macro in minwindef.h)
static inline bool approx(double a, double b, double eps = 1e-5)
{
    return (fabs(a - b) < eps);
}

static inline bool approx(const chai3d::cVector3d& a, const chai3d::cVector3d& b, double eps = 1e-5)
{
    return ((a - b).length() < eps);
}

static inline bool approx(const chai3d::cMatrix3d& a, const chai3d::cMatrix3d& b, double eps = 1e-5)
{
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            if (fabs(a(i,j) - b(i,j)) >= eps) return (false);
    return (true);
}

// Random draws are sequenced explicitly: the evaluation order of function
// arguments is unspecified, so cVector3d(d(rng), d(rng), d(rng)) could give
// different vectors with different compilers.
template <class D>
static inline chai3d::cVector3d randomVector(std::mt19937& a_rng, D& a_dist)
{
    const double x = a_dist(a_rng);
    const double y = a_dist(a_rng);
    const double z = a_dist(a_rng);
    return (chai3d::cVector3d(x, y, z));
}

// Uniform random point in the box [-a_half, a_half].
static inline chai3d::cVector3d randomInBox(std::mt19937& a_rng, const chai3d::cVector3d& a_half)
{
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    const chai3d::cVector3d r = randomVector(a_rng, uni);
    return (chai3d::cVector3d(r(0) * a_half(0), r(1) * a_half(1), r(2) * a_half(2)));
}

// Uniformly distributed random rotation.
static inline chai3d::cMatrix3d randomRotation(std::mt19937& a_rng)
{
    std::normal_distribution<double> gauss(0.0, 1.0);
    const double w = gauss(a_rng);
    const double x = gauss(a_rng);
    const double y = gauss(a_rng);
    const double z = gauss(a_rng);
    chai3d::cQuaternion q(w, x, y, z);
    q.normalize();
    chai3d::cMatrix3d rot;
    q.toRotMat(rot);
    return (rot);
}

static inline int testResult()
{
    printf(g_fails ? "%d check(s) FAILED\n" : "All checks passed\n", g_fails);
    return (g_fails ? 1 : 0);
}

//------------------------------------------------------------------------------
#endif
//------------------------------------------------------------------------------
