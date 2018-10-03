/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:17
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/licprocessor.h>
#include <queue>
#include <lablic/integrator.h>
#include <lablic/interpolator.h>
#include <inviwo/core/datastructures/volume/volumeram.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo LICProcessor::processorInfo_{
    "org.inviwo.LICProcessor",  // Class identifier
    "LICProcessor",             // Display name
    "KTH Labs",                 // Category
    CodeState::Experimental,    // Code state
    Tags::None,                 // Tags
};

const ProcessorInfo LICProcessor::getProcessorInfo() const { return processorInfo_; }

LICProcessor::LICProcessor()
    : Processor()
    , volumeIn_("volIn")
    , noiseTexIn_("noiseTexIn")
    , licOut_("licOut")
	, fast("fast", "Fast LIC", false)
	, kernelSize("kernelSize", "Kernel Size", 15, 3, 1000)
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
	addProperty(fast);
	addProperty(kernelSize);
}

void LICProcessor::process() {
    // Get input
    if (!volumeIn_.hasData()) {
        return;
    }

    if (!noiseTexIn_.hasData()) {
        return;
    }

    auto vol = volumeIn_.getData();
    vectorFieldDims_ = vol->getDimensions();
    auto vr = vol->getRepresentation<VolumeRAM>();
	auto dims = vr->getDimensions();

    // An accessible form of on image is retrieved analogous to a volume
    auto tex = noiseTexIn_.getData();
    texDims_ = tex->getDimensions();
    auto tr = tex->getRepresentation<ImageRAM>();

    // Prepare the output, it has the same dimensions and datatype as the output
    // and an editable version is retrieved analogous to a volume
    auto outImage = tex->clone();
    auto outLayer = outImage->getColorLayer();
    auto lr = outLayer->getEditableRepresentation<LayerRAM>();

    // To access the image at a floating point position, you can call
    //      Interpolator::sampleFromGrayscaleImage(tr, somePos)
	float pixelLength = (float) dims.x / texDims_.x;
	float pixelHeight = (float) dims.y / texDims_.y;
	float stepSize = pixelLength < pixelHeight ? pixelLength : pixelHeight;
    bool ** exploredPixels = new bool*[texDims_.x];
    for (unsigned int i = 0; i < texDims_.x; i++)
        exploredPixels[i] = new bool[texDims_.y];
    for (unsigned int j = 0; j < texDims_.y; j++) {
        for (unsigned int i = 0; i < texDims_.x; i++) {
			if (exploredPixels[i][j] == true) continue;
            // Create streamline
			auto vertices = Integrator::createStreamLine(vec2(i * stepSize, j * stepSize), vol.get(), 1000, stepSize);
            float currentSum = 0;
            std::queue<float> valueQueue;
            unsigned int queueHeader = 0;
            // Init first values in queue
            while(queueHeader < kernelSize.get() / 2) {
                float val = Interpolator::sampleFromGrayscaleImage(tr, vec2(vertices[queueHeader].x / stepSize, vertices[queueHeader].y / stepSize));
                currentSum += val;
                valueQueue.push(val);
                queueHeader++;
            }
            // find first point value
            float pointValue = currentSum / valueQueue.size();
            unsigned int pixelIdxX = vertices[0].x / stepSize;
            unsigned int pixelIdxY = vertices[0].y / stepSize;
            exploredPixels[pixelIdxX][pixelIdxY] = true;
            // Draw first point
            lr->setFromDVec4(size2_t(pixelIdxX, vpixelIdxY), dvec4(pointValue, pointValue, pointValue, 255));
			for(size_t l = 1; l < vertices.size(); l++) {
				vec2 vertex = vertices[l];
				pixelIdxX = vertex.x / stepSize;
				pixelIdxY = vertex.y / stepSize;
				exploredPixels[pixelIdxX][pixelIdxY] = true;
                if(valueQueue.size() > kernelSize.get()) {
                    currentSum -= valueQueue.front();
                    valueQueue.pop();
                }
                if(queueHeader < vertices.size()) {
                    float val = Interpolator::sampleFromGrayscaleImage(tr, vec2(vertices[queueHeader].x / stepSize, vertices[queueHeader].y / stepSize));
                    currentSum += val;
                    valueQueue.push(val);
                    queueHeader++;
                } else {
                    currentSum -= valueQueue.front();
                    valueQueue.pop();
                }
                pointValue = currentSum / valueQueue.size();
                lr->setFromDVec4(size2_t(pixelIdxX, pixelIdxY), dvec4(pointValue, pointValue, pointValue, 255));

			}
        }
    }

    licOut_.setData(outImage);
}

}  // namespace inviwo
