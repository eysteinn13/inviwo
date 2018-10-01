/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Tuesday, September 19, 2017 - 15:08:33
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <labstreamlines/streamlineintegrator.h>
#include <labstreamlines/integrator.h>
#include <inviwo/core/util/utilities.h>
#include <inviwo/core/interaction/events/mouseevent.h>
#include <math.h>
#include <time.h>
#include <vector>

namespace inviwo {


// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo StreamlineIntegrator::processorInfo_{
    "org.inviwo.StreamlineIntegrator",  // Class identifier
    "Streamline Integrator",            // Display name
    "KTH Lab",                          // Category
    CodeState::Experimental,            // Code state
    Tags::None,                         // Tags
};

const ProcessorInfo StreamlineIntegrator::getProcessorInfo() const { return processorInfo_; }

StreamlineIntegrator::StreamlineIntegrator()
	: Processor()
	, inData("volIn")
	, outMesh("meshOut")
	, propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
	, propSeedMode("seedMode", "Seeds")
	, multi_seed_method("seed_method", "Seed Method")
	, direction("direction", "Direction", 1, -1, 1, 2)
	, num_steps("num_steps", "Number of steps", 30, 1, 100, 1)
	, step_size("step_size", "Step size", 0.1, 0.0001, 10, 0.00001)
	, arc_lenght("arc_lenght", "Arc Length", 0.5, 0.001, 100, 0.001)
	, num_seeds("num_seeds", "Number of seeds", 1, 1, 100, 1)
	, num_seeds_x("num_seeds_x", "Number of seeds", 1, 1, 20, 1)
	, num_seeds_y("num_seeds_y", "Number of seeds", 1, 1, 20, 1)
	, stop_threshold("stop_threshold", "Stop Threshold", 0.0, 0.001, 5, 0.001)
	, direction_field("direction_field", "Use direction field", false)
    // TODO: Initialize additional properties
    // propertyName("propertyIdentifier", "Display Name of the Propery",
    // default value (optional), minimum value (optional), maximum value (optional), increment
    // (optional)); propertyIdentifier cannot have spaces
    , mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
                     MouseButton::Left, MouseState::Press | MouseState::Move) {
    // Register Ports
    addPort(inData);
    addPort(outMesh);

    // Register Properties
    propSeedMode.addOption("one", "Single Start Point", 0);
    propSeedMode.addOption("multiple", "Multiple Seeds", 1);
    addProperty(propSeedMode);
	multi_seed_method.addOption("random", "Random Seeds", 0);
	multi_seed_method.addOption("uniform", "Uniform Seeds", 1);
	multi_seed_method.addOption("uniform", "Multitute based", 2);
	addProperty(multi_seed_method);
	
    addProperty(propStartPoint);
    addProperty(mouseMoveStart);
	addProperty(arc_lenght);
	addProperty(stop_threshold);
	addProperty(num_seeds);
	addProperty(num_seeds_x);
	addProperty(num_seeds_y);

	addProperty(step_size);
	addProperty(num_steps);
	addProperty(direction);
	addProperty(direction_field);

    // TODO: Register additional properties
    // addProperty(propertyName);

    // You can hide and show properties for a single seed and hide properties for multiple seeds (TODO)
    propSeedMode.onChange([this]() {
        if (propSeedMode.get() == 0) {
            util::show(propStartPoint, mouseMoveStart);
			util::hide(num_seeds);
            // util::hide(...)
        } else {
			util::show(num_seeds);
            util::hide(propStartPoint, mouseMoveStart);
            // util::show(...)
        }
    });

}

void StreamlineIntegrator::eventMoveStart(Event* event) {
    // Handle mouse interaction only if we
    // are in the mode with a single point
    if (propSeedMode.get() == 1) return;
    auto mouseEvent = static_cast<MouseEvent*>(event);
    vec2 mousePos = mouseEvent->posNormalized();
    // Denormalize to volume dimensions
    mousePos.x *= dims.x - 1;
    mousePos.y *= dims.y - 1;
    // Update starting point
    propStartPoint.set(mousePos);
    event->markAsUsed();
}

void StreamlineIntegrator::process() {
    // Get input
//    if (!inData.hasData()) {
//        return;
//    }
//    auto vol = inData.getData();
//
//    // Retreive data in a form that we can access it
//    auto vr = vol->getRepresentation<VolumeRAM>();
//    dims = vol->getDimensions();
//    // The start point should be inside the volume (set maximum to the upper right corner)
//    propStartPoint.setMaxValue(vec2(dims.x - 1, dims.y - 1));
//
//    auto mesh = std::make_shared<BasicMesh>();
//    std::vector<BasicMesh::Vertex> vertices;
//
//    if (propSeedMode.get() == 0) {
//        auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
//        auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
//        
//        vec2 startPoint = propStartPoint.get();
//        vertices.push_back({vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
//                            vec3(0), vec3(0), vec4(0, 0, 0, 1)});
//        indexBufferPoints->add(static_cast<std::uint32_t>(0));        
//        createStreamLine(startPoint, vr, indexBufferRK, indexBufferPoints, vertices);
//        
//
//    } else {
//        if(multi_seed_method.get() == 0)
//        {
//            srand(time(NULL));
//            for (int i = 0; i < num_seeds.get(); i++)
//            {
//                auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
//                auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
//                // Draw start point
//                vec2 min_point = propStartPoint.getMinValue();
//                vec2 max_point = propStartPoint.getMaxValue();
//            
//                float random = ((float)rand()) / (float)RAND_MAX;
//                float diff = max_point.x - min_point.x;
//                float r = random * diff;
//                float rand_x = min_point.x + r;
//
//                random = ((float)rand()) / (float)RAND_MAX;
//                diff = max_point.y - min_point.y;
//                r = random * diff;
//                float rand_y = min_point.y + r;
//
//                vec2 startPoint(rand_x, rand_y);
//
//                indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
//                vertices.push_back({ vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
//                    vec3(0), vec3(0), vec4(0, 0, 0, 1) });            
//
//                createStreamLine(startPoint, vr, indexBufferRK, indexBufferPoints, vertices);
//            
//            }
//        }
//        else if (multi_seed_method.get() == 1)
//        {
//            vec2 min_point = propStartPoint.getMinValue();
//            vec2 max_point = propStartPoint.getMaxValue();
//
//            float increment_x = (max_point.x - min_point.x) / (num_seeds_x.get() + 1);
//            float increment_y = (max_point.y - min_point.y) / (num_seeds_y.get() + 1);
//
//            float x = min_point.x + increment_x;            
//            while(x < max_point.x)
//            //for (int x = min_point.x + increment_x; x < max_point.x; x += increment_x)
//            {
//                float y = min_point.y + increment_y;
//                while(y < max_point.y)
//                //for (int y = min_point.y + increment_y; y < max_point.y; y += increment_y)
//                {
//                    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
//                    auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
//                    vec2 startPoint(x, y);
//
//                    indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
//                    vertices.push_back({ vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
//                        vec3(0), vec3(0), vec4(0, 0, 0, 1) });
//
//                    createStreamLine(startPoint, vr, indexBufferRK, indexBufferPoints, vertices);
//                    y += increment_y;
//                }
//                x += increment_x;
//            }
//        }
//
//        else if (multi_seed_method.get() == 2)
//        {
//            srand(time(NULL));
//            int n = 15;
//            int sample_size = 20;
//            int points_sampled = 0;
//            vec2 min_point = propStartPoint.getMinValue();
//            vec2 max_point = propStartPoint.getMaxValue();
//            float increment_x = (max_point.x - min_point.x) / (sample_size + 1);
//            float increment_y = (max_point.y - min_point.y) / (sample_size + 1);
//            
//            float total_magnitude = 0;
//            int count_points = 0;
//            std::vector<float> magnitutes(676);
//            std::vector<vec2> vectors(676);
//            for (float x = increment_x; x < dims.x; x += increment_x)
//            {
//                for (float y = increment_y; y < dims.y; y += increment_y)
//                {
//                    vec2 sample = Integrator::sampleFromField(vr, dims, vec2(x, y), false);
//                    float size = sqrt(pow(sample[0], 2) + pow(sample[1], 2));
//                    vectors[count_points] = vec2(x, y);
//                    magnitutes[count_points++] = size;                            
//                    total_magnitude += size;
//                }
//            }
//
//            for (int i = 0; i < count_points; i++)
//            {
//                magnitutes[i] = (int) (magnitutes[i] * 1000 / total_magnitude);
//            }
//            std::vector<vec2> vec(n);
//            std::default_random_engine generator;
//            std::discrete_distribution<int> distribution(magnitutes.begin(), magnitutes.end());
//            std::vector<int> indices(vec.size());
//            std::generate(indices.begin(), indices.end(), [&generator, &distribution]() { return distribution(generator); });
//            std::transform(indices.begin(), indices.end(), vec.begin(), [&vectors](int index) { return vectors[index]; });
//            //for (int i = 0; i < indices.size(); i++)
//            for (int i = 0; i < vec.size(); i++)
//            {
//                auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);
//                auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
//                //vec2 startPoint(vectors[indices[i]]);
//                vec2 startPoint = vec[i];
//
//                indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
//                vertices.push_back({ vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
//                    vec3(0), vec3(0), vec4(0, 0, 0, 1) });
//
//                createStreamLine(startPoint, vr, indexBufferRK, indexBufferPoints, vertices);
//            }
//        }
//    }
//
//    mesh->addVertices(vertices);
//    outMesh.setData(mesh);
}

void StreamlineIntegrator::createStreamLine(const vec2 & startPoint, const VolumeRAM* vr, 
											IndexBufferRAM* indexBufferRK, IndexBufferRAM* indexBufferPoints,
											std::vector<BasicMesh::Vertex>& vertices)
{
	bool threshold = false;
	float arc_length_so_far = 0.0f;
	vec2 prev_point = startPoint;
	for (int i = 0; i < num_steps.get(); i++)
	{
		vec2 next_point = Integrator::RK4(vr, dims, prev_point, step_size.get(),
			direction.get(), direction_field.get(),
			threshold, stop_threshold.get());

		float x1 = next_point.x;
		float x2 = next_point[0];
		float y1 = next_point.y;
		float y2 = next_point[1];

		arc_length_so_far += (sqrt(pow(next_point.x - prev_point.x, 2) + pow(next_point.y - prev_point.y, 2)));
		bool out_of_bounds = (next_point[0] < 0 || next_point[0] > dims[0] - 1 || next_point[1] < 0 || next_point[1] > dims[1] - 1);

		if (threshold)
			LogProcessorInfo("threshold at " << next_point);

		if (arc_length_so_far >= arc_lenght || out_of_bounds || threshold)
			break;

		drawLineSegment(vec2(prev_point.x / (dims.x - 1), prev_point.y / (dims.y - 1)),
			vec2(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1)),
			vec4(1, 0, 0, 1), indexBufferRK, vertices);

		indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
		vertices.push_back({ vec3(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1), 0),
			vec3(0), vec3(0), vec4(1, 0, 0, 1) });

		prev_point = next_point;
	}
}

void StreamlineIntegrator::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
	IndexBufferRAM* indexBuffer,
	std::vector<BasicMesh::Vertex>& vertices) {
	// Add first vertex
	indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
	// A vertex has a position, a normal, a texture coordinate and a color
	// we do not use normal or texture coordinate, but still have to specify them
	vertices.push_back({ vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color });
	// Add second vertex
	indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
	vertices.push_back({ vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color });
}
}  // namespace inviwo
