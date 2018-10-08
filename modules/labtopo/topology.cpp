/*********************************************************************
*  Author  : Anke Friederici
*
*  Project : KTH Inviwo Modules
*
*  License : Follows the Inviwo BSD license model
**********************************************************************/

#include <inviwo/core/datastructures/geometry/basicmesh.h>
#include <inviwo/core/datastructures/volume/volumeram.h>
#include <labtopo/integrator.h>
#include <labtopo/interpolator.h>
#include <labtopo/topology.h>
#include <labtopo/utils/gradients.h>
#include <math.h>

namespace inviwo
{

const vec4 Topology::ColorsCP[6] =
    {
        vec4(1, 1, 0, 1),  // Saddle
        vec4(0, 0, 1, 1),  // AttractingNode
        vec4(1, 0, 0, 1),  // RepellingNode
        vec4(0.5, 0, 1, 1),// AttractingFocus
        vec4(1, 0.5, 0, 1),// RepellingFocus
        vec4(0, 1, 0, 1)   // Center
};

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo Topology::processorInfo_{
    "org.inviwo.Topology",  // Class identifier
    "Vector Field Topology",// Display name
    "KTH Lab",              // Category
    CodeState::Experimental,// Code state
    Tags::None,             // Tags
};

const ProcessorInfo Topology::getProcessorInfo() const
{
    return processorInfo_;
}

Topology::Topology()
    : Processor(), outMesh("meshOut"), inData("inData")
	, squareSizeThreshold("squareSizeThreshold", "Square Size Threshold", 0.3, 0.0001, 1, 0.0001)
// TODO: Initialize additional properties
// propertyName("propertyIdentifier", "Display Name of the Propery",
// default value (optional), minimum value (optional), maximum value (optional), increment (optional));
// propertyIdentifier cannot have spaces
{
    // Register Ports
    addPort(outMesh);
    addPort(inData);

    // TODO: Register additional properties
    // addProperty(propertyName);
	addProperty(squareSizeThreshold);
}

void Topology::process()
{
    // Get input
    if (!inData.hasData())
    {
        return;
    }
    auto vol = inData.getData();

    // Retreive data in a form that we can access it
    const VolumeRAM* vr = vol->getRepresentation<VolumeRAM>();
    uvec3 dims = vr->getDimensions();

    // Initialize mesh, vertices and index buffers for the two streamlines and the
    auto mesh = std::make_shared<BasicMesh>();
    std::vector<BasicMesh::Vertex> vertices;
    // Either add all line segments to this index buffer (one large buffer),
    // or use several index buffers with connectivity type adjacency.
    auto indexBufferSeparatrices = mesh->addIndexBuffer(DrawType::Lines, ConnectivityType::None);
    auto indexBufferPoints = mesh->addIndexBuffer(DrawType::Points, ConnectivityType::None);

    // TODO: Compute the topological skeleton of the input vector field.
    // Find the critical points and color them according to their type.
    // Integrate all separatrices.
    // You can use your previous integration code (copy it over or call it from <lablic/integrator.h>).

    // Looping through all values in the vector field.
	std::vector<vec2> criticalPoints;
	for (int y = 0; y < dims[1] - 1; ++y) {
		for (int x = 0; x < dims[0] - 1; ++x) {
			vec2 bottomLeft(x, y);
			vec2 bottomRight(x + 1, y);
			vec2 topLeft(x, y + 1);
			vec2 topRight(x + 1, y + 1);			
			zeroPossible(bottomLeft, bottomRight, topLeft, topRight, vol.get(), criticalPoints);
		}
	}

	LogProcessorInfo("Number of critical points: " << criticalPoints.size());
	
	for (auto point : criticalPoints)
	{
		indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));
		vertices.push_back({ vec3(point.x / (dims[0] - 1), point.y / (dims[1] - 1), 0), vec3(0), vec3(0), vec4(1, 0, 0, 1) });   
	}

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}

void Topology::zeroPossible(dvec2 bottomLeft, dvec2 bottomRight, dvec2 topLeft, dvec2 topRight, const Volume * vol, std::vector<vec2> & critPoints)
{
	float diffX = topRight.x - topLeft.x;
	float diffY = topRight.y - bottomRight.y;
	float squareSize = diffX * diffY;
	float offsetX = diffX / 2;
	float offsetY = diffY / 2;

	vec2 bottomLeftVal = Interpolator::sampleFromField(vol, bottomLeft);
	vec2 bottomRightVal = Interpolator::sampleFromField(vol, bottomRight);
	vec2 topLeftVal = Interpolator::sampleFromField(vol, topLeft);
	vec2 topRightVal = Interpolator::sampleFromField(vol, topRight);

	bool xPos = false;
	bool yPos = false;
	bool xNeg = false;
	bool yNeg = false;

	xPos = bottomLeftVal.x >= 0
		|| bottomRightVal.x >= 0
		|| topLeftVal.x >= 0
		|| topRightVal.x >= 0;

	yPos = bottomLeftVal.y >= 0
		|| bottomRightVal.y >= 0
		|| topLeftVal.y >= 0
		|| topRightVal.y >= 0;

	xNeg = bottomLeftVal.x < 0
		|| bottomRightVal.x < 0
		|| topLeftVal.x < 0
		|| topRightVal.x < 0;

	yNeg = bottomLeftVal.y < 0
		|| bottomRightVal.y < 0
		|| topLeftVal.y < 0
		|| topRightVal.y < 0;

	bool zeroInSquare = xPos && yPos && xNeg && yNeg;
	
	if (squareSize < squareSizeThreshold.get() && zeroInSquare)
	{
		dvec2 centerPoint = bottomLeft;
		centerPoint.x += offsetX;
		centerPoint.y += offsetY;

		critPoints.push_back(centerPoint);
		return;
	}
	else if(zeroInSquare)
	{
		dvec2 midPointLeft = bottomLeft;
		midPointLeft.y += offsetY;
		dvec2 midPointRight = bottomRight;
		midPointRight.y += offsetY;
		dvec2 midPointBottom = bottomLeft;
		midPointBottom.x += offsetX;
		dvec2 midPointTop = topLeft;
		midPointTop.x += offsetX;
		dvec2 centerPoint = bottomLeft;
		centerPoint.x += offsetX;
		centerPoint.y += offsetY;

		zeroPossible(bottomLeft, midPointBottom, midPointLeft, centerPoint, vol, critPoints);
		zeroPossible(midPointBottom, bottomRight, centerPoint, midPointRight, vol, critPoints);
		zeroPossible(midPointLeft, centerPoint, topLeft, midPointTop, vol, critPoints);
		zeroPossible(centerPoint, midPointRight, midPointTop, topRight, vol, critPoints);
	}
}

}// namespace
