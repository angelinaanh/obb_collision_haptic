//==============================================================================
/*
    OBB (Oriented Bounding Box) Tree extension for CHAI3D.

    \file       CCollisionOBBRefine.cpp
    \ingroup    obb
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "CCollisionOBBRefine.h"
//------------------------------------------------------------------------------
#include "math/CMaths.h"
//------------------------------------------------------------------------------
#include <cmath>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace chai3d {
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
namespace {
//------------------------------------------------------------------------------

//! Size and cost of the box around the points for a given frame.
struct cBoxCost
{
    double m_size[3];
    double m_volume;
    double m_area;
};


//! Projects the points on the axes and returns the size, volume and area of the box.
cBoxCost evaluate(const std::vector<cVector3d>& a_points, const cVector3d a_axes[3])
{
    // projections relative to the first point, which therefore projects to 0
    const cVector3d& ref = a_points[0];
    double lo[3] = { 0.0, 0.0, 0.0 };
    double hi[3] = { 0.0, 0.0, 0.0 };
    for (size_t i=1; i<a_points.size(); i++)
    {
        const cVector3d d = a_points[i] - ref;
        for (int k=0; k<3; k++)
        {
            const double t = cDot(d, a_axes[k]);
            if (t < lo[k]) lo[k] = t;
            if (t > hi[k]) hi[k] = t;
        }
    }

    cBoxCost cost;
    for (int k=0; k<3; k++) cost.m_size[k] = hi[k] - lo[k];
    cost.m_volume = cost.m_size[0] * cost.m_size[1] * cost.m_size[2];
    cost.m_area = 2.0 * (cost.m_size[0] * cost.m_size[1] +
                         cost.m_size[1] * cost.m_size[2] +
                         cost.m_size[2] * cost.m_size[0]);
    return (cost);
}


//! __true__ if a_candidate is a smaller box than a_current (volume first, then area).
bool isBetter(const cBoxCost& a_candidate, const cBoxCost& a_current,
              double a_volumeTolerance, double a_areaTolerance)
{
    if (a_candidate.m_volume < a_current.m_volume - a_volumeTolerance) return (true);
    if ((a_candidate.m_volume <= a_current.m_volume + a_volumeTolerance) &&
        (a_candidate.m_area < a_current.m_area - a_areaTolerance)) return (true);
    return (false);
}


//! Rotates the frame about its axis a_index; the frame stays orthonormal and right-handed.
void rotateAbout(cVector3d a_axes[3], int a_index, double a_angleRad)
{
    const int a = (a_index + 1) % 3;
    const int b = (a_index + 2) % 3;
    const double c = cos(a_angleRad);
    const double s = sin(a_angleRad);
    const cVector3d ua =  c * a_axes[a] + s * a_axes[b];
    const cVector3d ub = -s * a_axes[a] + c * a_axes[b];
    a_axes[a] = ua;
    a_axes[b] = ub;
}


//! Makes the frame exactly orthonormal and right-handed (Gram-Schmidt).
void orthonormalize(cVector3d a_axes[3])
{
    a_axes[0].normalize();
    a_axes[1] = cNormalize(a_axes[1] - cDot(a_axes[0], a_axes[1]) * a_axes[0]);
    a_axes[2] = cCross(a_axes[0], a_axes[1]);
}

//------------------------------------------------------------------------------
} // anonymous namespace
//------------------------------------------------------------------------------


//==============================================================================
/*!
    This function computes the size of the tightest box with the given axes
    around a point set.
*/
//==============================================================================
cVector3d cComputeBoxSize(const std::vector<cVector3d>& a_points,
                          const cMatrix3d& a_axes)
{
    if (a_points.empty())
    {
        return (cVector3d(0.0, 0.0, 0.0));
    }

    const cVector3d axes[3] = { a_axes.getCol0(), a_axes.getCol1(), a_axes.getCol2() };
    const cBoxCost cost = evaluate(a_points, axes);
    return (cVector3d(cost.m_size[0], cost.m_size[1], cost.m_size[2]));
}


//==============================================================================
/*!
    This function refines the axes of an OBB by greedy rotations about each of
    its axes, keeping every rotation that reduces the volume of the box.
*/
//==============================================================================
cMatrix3d cRefineOBBAxes(const std::vector<cVector3d>& a_points,
                         const cMatrix3d& a_axes,
                         const cOBBRefineSettings& a_settings,
                         cOBBRefineStats* a_stats)
{
    cOBBRefineStats stats;

    cVector3d u[3] = { a_axes.getCol0(), a_axes.getCol1(), a_axes.getCol2() };
    orthonormalize(u);

    if (a_points.empty())
    {
        if (a_stats != NULL) *a_stats = stats;
        return (cMatrix3d(u[0], u[1], u[2]));
    }

    cBoxCost current = evaluate(a_points, u);
    stats.m_initialVolume = current.m_volume;
    stats.m_evaluations = 1;

    // differences below these values are round-off, not improvements
    const double L = cMax(current.m_size[0], cMax(current.m_size[1], current.m_size[2]));
    const double volumeTolerance = 1e-12 * L * L * L;
    const double areaTolerance = 1e-12 * L * L;

    if (a_settings.m_enabled)
    {
        for (int iteration=0; iteration<a_settings.m_maxIterations; iteration++)
        {
            bool improved = false;

            for (int k=0; k<3; k++)
            {
                for (size_t s=0; s<a_settings.m_anglesDeg.size(); s++)
                {
                    const double angle = cDegToRad(a_settings.m_anglesDeg[s]);

                    // try +angle and -angle, keep the better one if it improves the box
                    cBoxCost best = current;
                    cVector3d bestAxes[3];
                    bool found = false;
                    for (int sign=-1; sign<=1; sign+=2)
                    {
                        cVector3d candidate[3] = { u[0], u[1], u[2] };
                        rotateAbout(candidate, k, sign * angle);
                        const cBoxCost cost = evaluate(a_points, candidate);
                        stats.m_evaluations++;

                        if (isBetter(cost, best, volumeTolerance, areaTolerance))
                        {
                            best = cost;
                            for (int j=0; j<3; j++) bestAxes[j] = candidate[j];
                            found = true;
                        }
                    }

                    if (found)
                    {
                        for (int j=0; j<3; j++) u[j] = bestAxes[j];
                        current = best;
                        improved = true;
                        stats.m_acceptedRotations++;
                    }
                }
            }

            stats.m_iterations++;
            if (!improved) break;
        }
    }

    orthonormalize(u);
    stats.m_finalVolume = current.m_volume;
    if (a_stats != NULL) *a_stats = stats;

    return (cMatrix3d(u[0], u[1], u[2]));
}

//------------------------------------------------------------------------------
} // namespace chai3d
//------------------------------------------------------------------------------
