/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Wednesday, September 20, 2017 - 12:04:15
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/integrator.h>
#include <lablic/interpolator.h>

namespace inviwo {

Integrator::Integrator() {}
    
vec2 Integrator::RK4(const Volume* vol, const vec2& position, float step_size)
{
    vec2 v1 = Interpolator::sampleFromField(vol, position);
    vec2 v2 = Interpolator::sampleFromField(vol, position + (step_size / 2) * v1);
    vec2 v3 = Interpolator::sampleFromField(vol, position + (step_size / 2) * v2);
    vec2 v4 = Interpolator::sampleFromField(vol, position + (step_size * v3));
    
    float step_x = v1.x / 6 + v2.x / 3 + v3.x / 3 + v4.x / 6;
    float step_y = v1.y / 6 + v2.y / 3 + v3.y / 3 + v4.y / 6;
    vec2 next_point(position.x + (step_size * step_x),
                    position.y + (step_size * step_y));
    return next_point;
}
    
std::vector<vec2> Integrator::createStreamLine(const vec2 & startPoint, const Volume* vol, float arcLength, float stepSize)
{
    auto vr = vol->getRepresentation<VolumeRAM>();
    auto dims = vr->getDimensions();
    bool outOfBounds = false;
    vec2 prevPoint, nextPoint = startPoint;
    std::vector<vec2> vertices;
    for (int i = 0; i < 0.5 * arcLength; i += stepSize)
    {
        nextPoint = RK4(vol, prevPoint, stepSize);
        outOfBounds = (nextPoint[0] < 0 || nextPoint[0] > dims[0] - 1 || nextPoint[1] < 0 || nextPoint[1] > dims[1] - 1);

        if (outOfBounds)
            break;
        vertices.push_back(nextPoint);
        prevPoint = nextPoint;
    }
    prevPoint = startPoint;
    for (int i = 0; i < 0.5 * arcLength; i += stepSize, prevPoint = nextPoint)
    {
        nextPoint = RK4(vol, prevPoint, -stepSize);
        outOfBounds = (nextPoint[0] < 0 || nextPoint[0] > dims[0] - 1 || nextPoint[1] < 0 || nextPoint[1] > dims[1] - 1);
        
        if (outOfBounds)
            break;
        vertices.push_back(nextPoint);
    }
    return vertices;
}


// TODO: Implementation of the functions from last lab
// HINT: There is a change in sampleFromField():
//       Interpolator::sampleFromField(vol.get(), somePosition);

}  // namespace inviwo
