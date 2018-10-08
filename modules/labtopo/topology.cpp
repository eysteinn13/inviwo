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
#include <utility>

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
	, squareSizeThreshold("squareSizeThreshold", "Square Size Threshold", 0.3, 0.00000001, 1, 0.00000001)
	, zeroTolerance("zeroTolerance", "Zero Tolerance", 1e-4, 1e-7, 1e-1, 1e-1)
	, stepSize("stepSize", "Step Size", 1, 0.00001, 5, 0.00001)
	, threshold("threshold", "Threshold", 1, 0.00001, 5, 0.00001)
	, numSteps("numSteps", "Number of Steps", 100, 1, 500, 1)
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
	addProperty(zeroTolerance);
	addProperty(stepSize);
	addProperty(threshold);
	addProperty(numSteps);
}

void Topology::process()
{
    // Get input
    if (!inData.hasData()) {
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
	for (unsigned int y = 0; y < dims[1] - 1; ++y) {
		for (unsigned int x = 0; x < dims[0] - 1; ++x) {
			vec2 bottomLeft(x, y);
			vec2 bottomRight(x + 1, y);
			vec2 topLeft(x, y + 1);
			vec2 topRight(x + 1, y + 1);			
			findCriticalPoints(bottomLeft, bottomRight, topLeft, topRight, vol.get(), criticalPoints);
		}
	}

	LogProcessorInfo("Number of critical points: " << criticalPoints.size());
	
	for (auto point : criticalPoints) {
		vec4 color = getCritPointColor(point, vol.get(), indexBufferSeparatrices, vertices);
		indexBufferPoints->add(static_cast<std::uint32_t>(vertices.size()));		
		vertices.push_back({	vec3(point.x / (dims[0] - 1), point.y / (dims[1] - 1), 0), 
								vec3(0), vec3(0),  color});
	}

    mesh->addVertices(vertices);
    outMesh.setData(mesh);
}


vec4 Topology::getCritPointColor(vec2 point, const Volume * vol, IndexBufferRAM* indexBufferRK, 
								std::vector<BasicMesh::Vertex>& vertices)
{
	//StreamlineIntegrator * integrator = new StreamlineIntegrator();
	mat2 jacobian = Interpolator::sampleJacobian(vol, point);
	auto eigenRes = util::eigenAnalysis(jacobian);
	vec2 imaginaries = eigenRes.eigenvaluesIm;
	vec2 reals = eigenRes.eigenvaluesRe;

	if (approxEq(imaginaries[0], 0) && approxEq(imaginaries[1], 0)) {
		if ((reals[0] < 0 && reals[1] > 0) || (reals[1] < 0 && reals[0] > 0))
		{ 			
			vec2 startPoint1 = point + 1e-2f * eigenRes.eigenvectors[0];
			vec2 startPoint2 = point - 1e-2f  * eigenRes.eigenvectors[0];
			vec2 startPoint3 = point + 1e-2f  * eigenRes.eigenvectors[1];
			vec2 startPoint4 = point - 1e-2f * eigenRes.eigenvectors[1];

			int orientation1 = eigenRes.eigenvaluesRe[0] / abs(eigenRes.eigenvaluesRe[0]);
			int orientation2 = eigenRes.eigenvaluesRe[1] / abs(eigenRes.eigenvaluesRe[1]);

			createStreamLine(startPoint1, vol, orientation1, indexBufferRK, vertices);			
			createStreamLine(startPoint3, vol, orientation1, indexBufferRK, vertices);

			createStreamLine(startPoint2, vol, orientation2, indexBufferRK, vertices);
			createStreamLine(startPoint4, vol, orientation2, indexBufferRK, vertices);
			return ColorsCP[0];
		}
		if (reals[0] < 0 && reals[1] < 0)
			return ColorsCP[1];
		if (reals[0] > 0 && reals[1] > 0)
			return ColorsCP[2];
	}
	else if (approxEq(imaginaries[0], -1 * imaginaries[1]) && !approxEq(imaginaries[0], 0)) {
		if (approxEq(reals[0], reals[1]) && reals[0] < 0)
			return ColorsCP[3];
        if (approxEq(reals[0], 0) && approxEq(reals[1], 0))
            return ColorsCP[5];
		if (reals[0] > 0 && reals[1] > 0)
			return ColorsCP[4];	}
	return vec4(0, 0, 0, 1);
}


void Topology::createStreamLine(const vec2 & startPoint, const Volume * vol, int direction,
								IndexBufferRAM* indexBufferRK, std::vector<BasicMesh::Vertex>& vertices)
{
	bool reachedThreshold = false;
	auto dims = vol->getDimensions();
	
	vec2 prev_point = startPoint;
	for (int i = 0; i < numSteps.get(); i++)
	{
		vec2 next_point = Integrator::RK4(vol, prev_point, stepSize.get(), direction, reachedThreshold, threshold.get());

		bool out_of_bounds = (next_point[0] < 0 || next_point[0] > dims[0] - 1 || next_point[1] < 0 || next_point[1] > dims[1] - 1);
		if (out_of_bounds || reachedThreshold)
			break;

		drawLineSegment(vec2(prev_point.x / (dims.x - 1), prev_point.y / (dims.y - 1)),
						vec2(next_point.x / (dims.x - 1), next_point.y / (dims.y - 1)),
						vec4(1, 0, 0, 1), indexBufferRK, vertices);
		prev_point = next_point;
	}
}


void Topology::drawLineSegment(const vec2& v1, const vec2& v2, const vec4& color,
	IndexBufferRAM* indexBuffer, std::vector<BasicMesh::Vertex>& vertices) {
	// Add first vertex
	indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
	// A vertex has a position, a normal, a texture coordinate and a color
	// we do not use normal or texture coordinate, but still have to specify them
	vertices.push_back({ vec3(v1[0], v1[1], 0), vec3(0, 0, 1), vec3(v1[0], v1[1], 0), color });
	// Add second vertex
	indexBuffer->add(static_cast<std::uint32_t>(vertices.size()));
	vertices.push_back({ vec3(v2[0], v2[1], 0), vec3(0, 0, 1), vec3(v2[0], v2[1], 0), color });
}

bool Topology::approxEq(float a, float b)
{
	return a > b - zeroTolerance.get() && a < b + zeroTolerance.get();
}

void Topology::findCriticalPoints(vec2 bottomLeft, vec2 bottomRight, vec2 topLeft, vec2 topRight, const Volume * vol, std::vector<vec2> & critPoints)
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

	bool xPos, yPos, xNeg, yNeg = false;
	xPos = bottomLeftVal.x > 0
		|| bottomRightVal.x > 0
		|| topLeftVal.x > 0
		|| topRightVal.x > 0;

	yPos = bottomLeftVal.y > 0
		|| bottomRightVal.y > 0
		|| topLeftVal.y > 0
		|| topRightVal.y > 0;

	xNeg = bottomLeftVal.x < 0
		|| bottomRightVal.x < 0
		|| topLeftVal.x < 0
		|| topRightVal.x < 0;

	yNeg = bottomLeftVal.y < 0
		|| bottomRightVal.y < 0
		|| topLeftVal.y < 0
		|| topRightVal.y < 0;

	bool zeroInSquare = xPos && yPos && xNeg && yNeg;
	
	if (squareSize < squareSizeThreshold.get() && zeroInSquare) {
		dvec2 centerPoint = bottomLeft;
		centerPoint.x += offsetX;
		centerPoint.y += offsetY;

		critPoints.push_back(centerPoint);
		return;
	}
	else if(zeroInSquare) {
		vec2 midPointLeft = bottomLeft;
		midPointLeft.y += offsetY;
		vec2 midPointRight = bottomRight;
		midPointRight.y += offsetY;
		vec2 midPointBottom = bottomLeft;
		midPointBottom.x += offsetX;
		vec2 midPointTop = topLeft;
		midPointTop.x += offsetX;
		vec2 centerPoint = bottomLeft;
		centerPoint.x += offsetX;
		centerPoint.y += offsetY;

		findCriticalPoints(bottomLeft, midPointBottom, midPointLeft, centerPoint, vol, critPoints);
		findCriticalPoints(midPointBottom, bottomRight, centerPoint, midPointRight, vol, critPoints);
		findCriticalPoints(midPointLeft, centerPoint, topLeft, midPointTop, vol, critPoints);
		findCriticalPoints(centerPoint, midPointRight, midPointTop, topRight, vol, critPoints);
	}
}

}// namespace
