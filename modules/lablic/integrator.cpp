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
    
vec2 Integrator::RK4(const Volume* vol, const vec2& position, float step_size,
                     int direction, bool dir_field, bool&threshold,
                     const float& threshold_value)
{
    vec2 v1 = Interpolator::sampleFromField(vol, position);
    vec2 v2 = Interpolator::sampleFromField(vol, position + (step_size / 2) * v1);
    vec2 v3 = Interpolator::sampleFromField(vol, position + (step_size / 2) * v2);
    vec2 v4 = Interpolator::sampleFromField(vol, position + (step_size * v3));
    
    float step_x = direction * (v1.x / 6 + v2.x / 3 + v3.x / 3 + v4.x / 6);
    float step_y = direction * (v1.y / 6 + v2.y / 3 + v3.y / 3 + v4.y / 6);
    float denom = sqrt(pow(step_x, 2) + pow(step_y, 2));
    
    if (denom < threshold_value || denom == 0)
        threshold = true;
    
    vec2 next_point(position.x + (step_size * step_x),
                    position.y + (step_size * step_y));
    return next_point;
}

// TODO: Implementation of the functions from last lab
// HINT: There is a change in sampleFromField():
//       Interpolator::sampleFromField(vol.get(), somePosition);

}  // namespace inviwo
