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
    
vec2 Integrator::RK4(const Volume* vol, const vec2& position, float step_size, bool & belowThreshold)
{   
	float magnitute = 0;
	vec2 v1 = Interpolator::sampleFromField(vol, position, belowThreshold, magnitute);
    vec2 v2 = Interpolator::sampleFromField(vol, position + (step_size / 2) * v1, belowThreshold, magnitute);
    vec2 v3 = Interpolator::sampleFromField(vol, position + (step_size / 2) * v2, belowThreshold, magnitute);
    vec2 v4 = Interpolator::sampleFromField(vol, position + (step_size * v3), belowThreshold, magnitute);
    
    float step_x = v1.x / 6 + v2.x / 3 + v3.x / 3 + v4.x / 6;
    float step_y = v1.y / 6 + v2.y / 3 + v3.y / 3 + v4.y / 6;
    vec2 next_point(position.x + (step_size * step_x),
                    position.y + (step_size * step_y));
    return next_point;
}

std::vector<vec2> Integrator::createStreamLineSlow(const vec2 & startPoint, const Volume* vol, float arcLength, float stepSize, float pixleL, float pixleH, size2_t dimsTex, float & largest)
{
	auto vr = vol->getRepresentation<VolumeRAM>();
	auto dims = vr->getDimensions();
	bool outOfBounds = false;
	bool belowThreshold = false;
	vec2 prevPoint, nextPoint = startPoint;
	std::vector<vec2> vertices;
	for (float i = 0; i < 0.5 * arcLength; i += stepSize)
	{
		nextPoint = RK4(vol, prevPoint, stepSize, belowThreshold);
		outOfBounds = (nextPoint[0] < 0 || nextPoint[0] > dims[0] - 1 || nextPoint[1] < 0 || nextPoint[1] > dims[1] - 1);
		bool pixelOffGrid = (((int)(nextPoint.x / pixleL) > dimsTex.x - 2) || ((int)(nextPoint.y / pixleH) > dimsTex.y - 2));
		if (outOfBounds || pixelOffGrid)
			break;

		bool threshold;
		float magnitutue = 0.0f;
		Interpolator::sampleFromField(vol, nextPoint, threshold, magnitutue);
		if (magnitutue > largest)
			largest = magnitutue;

		vertices.push_back(nextPoint);
		prevPoint = nextPoint;
	}
	prevPoint = startPoint;
	for (float i = 0; i < 0.5 * arcLength; i += stepSize)
	{
		nextPoint = RK4(vol, prevPoint, -stepSize, belowThreshold);
		outOfBounds = (nextPoint[0] < 0 || nextPoint[0] > dims[0] - 1 || nextPoint[1] < 0 || nextPoint[1] > dims[1] - 1);
		bool pixelOffGrid = (((int)(nextPoint.x / pixleL) > dimsTex.x - 2) || ((int)(nextPoint.y / pixleH) > dimsTex.y - 2));
		if (outOfBounds || pixelOffGrid)
			break;

		bool threshold;
		float magnitutue = 0.0f;
		Interpolator::sampleFromField(vol, nextPoint, threshold, magnitutue);
		if (magnitutue > largest)
			largest = magnitutue;

		vertices.push_back(nextPoint);
		prevPoint = nextPoint;
	}
	return vertices;
}

std::vector<vec2> Integrator::createStreamLine(const vec2 & startPoint, const Volume* vol, float arcLength, float stepSize, bool ** explored , float pixleL, float pixleH, float & largest)
{
    auto vr = vol->getRepresentation<VolumeRAM>();
    auto dims = vr->getDimensions();
    bool outOfBounds = false;
	bool belowThreshold = false;
    vec2 prevPoint, nextPoint = startPoint;
    int countSamePixel = 0;
    std::vector<vec2> verticesFront;
    verticesFront.push_back(startPoint);
    // TODO add threshold check
    for (float i = 0; i < 0.5 * arcLength; i += stepSize)
    {
        nextPoint = RK4(vol, prevPoint, stepSize, belowThreshold);
        outOfBounds = (nextPoint[0] < 0 || nextPoint[0] > dims[0] - 1 || nextPoint[1] < 0 || nextPoint[1] > dims[1] - 1);

        if (outOfBounds)
            break;
        prevPoint = nextPoint;
        if(!explored[(int)(nextPoint.x / pixleL)][(int)(nextPoint.y / pixleH)]) {
            explored[(int)(nextPoint.x / pixleL)][(int)(nextPoint.y / pixleH)] = true;
            countSamePixel = 0;
        } else {
            countSamePixel++;
            if(countSamePixel > 100) {
                break;
            }
        }
        bool threshold;
        float magnitutue = 0.0f;
        Interpolator::sampleFromField(vol, nextPoint, threshold, magnitutue);
        if (magnitutue > largest)
            largest = magnitutue;
        verticesFront.push_back(nextPoint);
    }
    prevPoint = startPoint;

	std::vector<vec2> vertices_back;
    countSamePixel = 0;
    for (float i = 0; i < 0.5 * arcLength; i += stepSize)
    {
        nextPoint = RK4(vol, prevPoint, -stepSize, belowThreshold);
        outOfBounds = (nextPoint[0] < 0 || nextPoint[0] > dims[0] - 1 || nextPoint[1] < 0 || nextPoint[1] > dims[1] - 1);
        if (outOfBounds)
            break;
        prevPoint = nextPoint;
        if(!explored[(int)(nextPoint.x / pixleL)][(int)(nextPoint.y / pixleH)]) {
            explored[(int)(nextPoint.x / pixleL)][(int)(nextPoint.y / pixleH)] = true;
            countSamePixel = 0;
        } else {
            countSamePixel++;
            if(countSamePixel > 100) {
                break;
            }
        }
        bool threshold;
        float magnitutue = 0.0f;
        Interpolator::sampleFromField(vol, nextPoint, threshold, magnitutue);
        if (magnitutue > largest)
            largest = magnitutue;
		vertices_back.push_back(nextPoint);
    }
	std::reverse(vertices_back.begin(), vertices_back.end());
	vertices_back.insert(vertices_back.end(), verticesFront.begin(), verticesFront.end());
    return vertices_back;
}


// TODO: Implementation of the functions from last lab
// HINT: There is a change in sampleFromField():
//       Interpolator::sampleFromField(vol.get(), somePosition);

}  // namespace inviwo
