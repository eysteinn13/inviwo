/*********************************************************************
*  Author  : Himangshu Saikia
*  Init    : Wednesday, September 20, 2017 - 12:04:15
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
*********************************************************************
*/

#include <labstreamlines/integrator.h>
#include <math.h>

namespace inviwo
{

	Integrator::Integrator()
	{
	}

	vec2 Integrator::sampleFromField(const VolumeRAM* vr, size3_t dims, const vec2& position, const bool& normalize)
	{
		// Sampled outside the domain!
		if (position[0] < 0 || position[0] > dims[0] - 1 || position[1] < 0 || position[1] > dims[1] - 1)
		{
			return vec2(0, 0);
		}

		int x0 = int(position[0]);
		int y0 = int(position[1]);

		// Leads to accessing only inside the volume
		// Coefficients computation takes care of using the correct values
		if (x0 == dims[0] - 1)
		{
			x0--;
		}
		if (y0 == dims[1] - 1)
		{
			y0--;
		}

		auto f00 = vr->getAsDVec2(size3_t(x0, y0, 0));
		auto f10 = vr->getAsDVec2(size3_t(x0 + 1, y0, 0));
		auto f01 = vr->getAsDVec2(size3_t(x0, y0 + 1, 0));
		auto f11 = vr->getAsDVec2(size3_t(x0 + 1, y0 + 1, 0));

		float x = position[0] - x0;
		float y = position[1] - y0;

		vec2 f;

		for (int i = 0; i < 2; i++)
		{
			f[i] = f00[i] * (1 - x) * (1 - y) + f01[i] * (1 - x) * y + f10[i] * x * (1 - y) + f11[i] * x * y;
		}
		if (normalize)
		{
			float denom = sqrt(pow(f[0], 2) + pow(f[1], 2));
			if (denom > 1e-7)
			{
				f[0] = f[0] / denom;
				f[1] = f[1] / denom;
			}
		}
		return f;
	}


	// TODO: Implement a single integration step here

	vec2 Integrator::Euler(const VolumeRAM* vr, size3_t dims, const vec2& position, float step_size, int direction)
	{
		vec2 sample = sampleFromField(vr, dims, position, direction);
		vec2 next_point(position.x + direction * (step_size * sample.x), 
						position.y + direction * (step_size * sample.y));
		return next_point;
	}

	vec2 Integrator::RK4(	const VolumeRAM* vr, size3_t dims, 
							const vec2& position, float step_size, 
							int direction, bool dir_field, bool&threshold,
							const float& threshold_value)
	{
		vec2 v1 = sampleFromField(vr, dims, position, dir_field);
		vec2 v2 = sampleFromField(vr, dims, position + (step_size / 2) * v1, dir_field);
		vec2 v3 = sampleFromField(vr, dims, position + (step_size / 2) * v2, dir_field);
		vec2 v4 = sampleFromField(vr, dims, position + (step_size * v3), dir_field);

		float step_x = direction * (v1.x / 6 + v2.x / 3 + v3.x / 3 + v4.x / 6);
		float step_y = direction * (v1.y / 6 + v2.y / 3 + v3.y / 3 + v4.y / 6);
		float denom = sqrt(pow(step_x, 2) + pow(step_y, 2));

		if (denom < threshold_value || denom == 0)
			threshold = true;

		vec2 next_point(position.x + (step_size * step_x),
						position.y + (step_size * step_y));
		return next_point;
	}


} // namespace

