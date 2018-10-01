/*********************************************************************
*  Author  : Himangshu Saikia, Wiebke Koepp
*  Init    : Tuesday, September 19, 2017 - 15:08:24
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
*********************************************************************
*/

#include <labstreamlines/eulerrk4comparison.h>
#include <labstreamlines/integrator.h>
#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/interaction/events/mouseevent.h>

namespace inviwo
{
	// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
	const ProcessorInfo EulerRK4Comparison::processorInfo_
	{
		"org.inviwo.EulerRK4Comparison",      // Class identifier
		"Euler RK4 Comparison",               // Display name
		"KTH Lab",                        // Category
		CodeState::Experimental,              // Code state
		Tags::None,                           // Tags
	};

	const ProcessorInfo EulerRK4Comparison::getProcessorInfo() const
	{
		return processorInfo_;
	}

	EulerRK4Comparison::EulerRK4Comparison()
		:Processor()
		, outMesh("meshOut")
		, inData("inData")
		, propStartPoint("startPoint", "Start Point", vec2(0.5, 0.5), vec2(0), vec2(1024), vec2(0.5))
		, num_steps_euler("num_steps_euler", "Number of steps Euler", 10, 100, 500, 1)
		, num_steps_kutta("num_steps_kutta", "Number of steps Kutta", 10, 30, 100, 1)
		, step_size_euler("step_size_euler", "Step size Euler", 0.1, 0.0001, 10, 0.00001)
		, step_size_kutta("step_size_kutta", "Step size Kutta", 0.1, 0.0001, 10, 0.00001)
		// TODO: Initialize additional properties
		// propertyName("propertyIdentifier", "Display Name of the Propery", 
		// default value (optional), minimum value (optional), maximum value (optional), increment (optional));
		// propertyIdentifier cannot have spaces
		, mouseMoveStart("mouseMoveStart", "Move Start", [this](Event* e) { eventMoveStart(e); },
			MouseButton::Left, MouseState::Press | MouseState::Move)
	{
		// Register Ports
		addPort(outMesh);
		addPort(inData);

		// Register Properties
		addProperty(propStartPoint);
		addProperty(mouseMoveStart);
		addProperty(num_steps_euler);
		addProperty(num_steps_kutta);
		addProperty(step_size_euler);
		addProperty(step_size_kutta);

		// TODO: Register additional properties
		// addProperty(propertyName);

	}

	void EulerRK4Comparison::eventMoveStart(Event* event)
	{
		auto mouseEvent = static_cast<MouseEvent*>(event);
		vec2 mousePos = mouseEvent->posNormalized();
		// Denormalize to volume dimensions
		mousePos.x *= dims.x - 1;
		mousePos.y *= dims.y - 1;
		// Update starting point
		propStartPoint.set(mousePos);
		event->markAsUsed();
	}

	void EulerRK4Comparison::process()
	{
		// Get input
		if (!inData.hasData())
		{
			return;
		}
		auto vol = inData.getData();

		// Retreive data in a form that we can access it
		const VolumeRAM* vr = vol->getRepresentation< VolumeRAM >();
		dims = vol->getDimensions();
		// The start point should be inside the volume (set maximum to the upper right corner)
		propStartPoint.setMaxValue(vec2(dims.x - 1, dims.y - 1));

		// Initialize mesh, vertices and index buffers for the two streamlines and the 
		auto mesh = std::make_shared<BasicMesh>();
		std::vector<BasicMesh::Vertex> vertices;
		auto indexBufferEuler = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		auto indexBufferRK = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::Strip);
		auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

		// Draw start point
		vec2 startPoint = propStartPoint.get();
		vertices.push_back({ vec3(startPoint.x / (dims.x - 1), startPoint.y / (dims.y - 1), 0),
			vec3(0), vec3(0), vec4(0, 0, 0, 1) });
		indexBufferPoints->add(static_cast<std::uint32_t>(0));
		indexBufferEuler->add(static_cast<std::uint32_t>(0));
		indexBufferRK->add(static_cast<std::uint32_t>(0));

		// TODO: Implement the Euler and Runge-Kutta of 4th order integration schemes
		// and then integrate forward for a specified number of integration steps and a given stepsize 
		// (these should be additional properties of the processor)
		LogProcessorInfo("Step size Euler" << step_size_euler.get());
		LogProcessorInfo("Step size Kutta" << step_size_kutta.get());
		vec2 prev_point = startPoint;

		for (int i = 0; i < num_steps_euler.get(); i++)
		{
			vec2 next_point = Integrator::Euler(vr, dims, prev_point, step_size_euler.get(), 1);
			drawLineSegment(vec2(prev_point.x / (dims.x - 1), prev_point.y / (dims.y - 1)),
				vec2(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1)),
				vec4(0, 0, 1, 1), indexBufferEuler, vertices);

			indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
			vertices.push_back({vec3(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1), 0),
								vec3(0), vec3(0), vec4(0, 0, 1, 1) });

			prev_point = next_point;
		}

		prev_point = startPoint;
		for (int i = 0; i < num_steps_kutta.get(); i++)
		{
			//indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
			//vertices.push_back({ vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color });

			bool temp = false;
			vec2 next_point = Integrator::RK4(vr, dims, prev_point, step_size_kutta.get(), 1, false, temp, 2.0);
			drawLineSegment(vec2(prev_point.x / (dims.x - 1), prev_point.y / (dims.y - 1)),
				vec2(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1)),
				vec4(1, 0, 0, 1), indexBufferRK, vertices);

			indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
			vertices.push_back({ vec3(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1), 0),
				vec3(0), vec3(0), vec4(1, 0, 0, 1) });

			prev_point = next_point;
		}

		// Integrator::Euler(vr, dims, startPoint, ...);
		// Integrator::Rk4(vr, dims, startPoint, ...);

		mesh->addVertices(vertices);
		outMesh.setData(mesh);
	}

	void EulerRK4Comparison::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
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

} // namespace

