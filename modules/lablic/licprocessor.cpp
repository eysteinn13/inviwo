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
	, kernelSize("kernelSize", "Kernel Size", 15, 3, 100)

// TODO: Register additional properties

{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
	addProperty(fast);
	addProperty(kernelSize);

    // TODO: Register additional properties

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

    // TODO: Implement LIC and FastLIC
    // This code instead sets all pixels to the same gray value
    // std::vector<std::vector<double>> licTexture(texDims_.x, std::vector<double>(texDims_.y, 127.0));

	float pixelLength = (float) dims.x / texDims_.x;
	float pixelHeight = (float) dims.y / texDims_.y;
	float stepSize = pixelLength < pixelHeight ? pixelLength : pixelHeight;
	bool ** exploredPixels = new bool*[texDims_.x];
	for (int i = 0; i < texDims_.x; i++)
	{
		exploredPixels[i] = new bool[texDims_.y];
		for (int j = 0; j < texDims_.y; j++)
			exploredPixels[i][j] = false;
	}

	float arcLength = kernelSize.get();
    for (auto j = 0; j < texDims_.y; j++) {
        for (auto i = 0; i < texDims_.x; i++) {

			if (exploredPixels[i][j] == true) continue;

			auto vertices = Integrator::createStreamLine(vec2(i * stepSize, j * stepSize), vol.get(), 1000, stepSize);						
			for(int i = 0; i < vertices.size(); i++)
			{
				vec2 vertex = vertices[i];
				int pixelIdxX = vertex.x / stepSize;
				int pixelIdxY = vertex.y / stepSize;
				exploredPixels[pixelIdxX][pixelIdxY] = true;

				int idx = i;
				int counter = 0;
				int stepCount = 0;
				float val = 0.0f;

				for (float k = 0; k < kernelSize.get() / 2; k += stepSize)
				{
					if (i + stepCount > vertices.size()) break;
					pixelIdxX = vertices[i + stepCount].x / stepSize;
					pixelIdxY = vertices[i + stepCount].y / stepSize;
					if (pixelIdxX < 0 || pixelIdxX > texDims_.x) break;
					if (pixelIdxY < 0 || pixelIdxY > texDims_.y) break;

					val += Interpolator::sampleFromGrayscaleImage(tr, vec2(pixelIdxX, pixelIdxY));
					counter++;
				}
				stepCount = 1;
				for (float k = stepSize; k < (kernelSize.get() - 1) / 2; k += stepSize)
				{
					if (i - stepCount < 0) break;
					pixelIdxX = vertices[i - stepCount].x / stepSize;
					pixelIdxY = vertices[i - stepCount].y / stepSize;
					if (pixelIdxX < 0 || pixelIdxX > texDims_.x) break;
					if (pixelIdxY < 0 || pixelIdxY > texDims_.y) break;

					val += Interpolator::sampleFromGrayscaleImage(tr, vec2(pixelIdxX, pixelIdxY));
					counter++;
				}
				val /= counter;
				lr->setFromDVec4(size2_t(i, j), dvec4(val, val, val, 255));			
			}
        }
    }

    licOut_.setData(outImage);
}

}  // namespace inviwo
