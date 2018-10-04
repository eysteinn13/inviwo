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
    , desiredU("desiredU", "U for contrast", 128, 0, 255)
    , desiredSigma("desiredSigma", "Sigma for contrast", 50, 0, 100, 0.001)
    , useContrast("useContrast", "Use contrast", false)
{
    // Register ports
    addPort(volumeIn_);
    addPort(noiseTexIn_);
    addPort(licOut_);

    // Register properties
	addProperty(fast);
	addProperty(kernelSize);
    addProperty(desiredU);
    addProperty(desiredSigma);
    addProperty(useContrast);
}

    void LICProcessor::contrast(LayerRAM * lr, const ImageRAM * tr){
        float u = 0;
        float P = 0;
        int n = 0;
        for (unsigned int j = 0; j < texDims_.y; j++) {
            for (unsigned int i = 0; i < texDims_.x; i++) {
                float val = Interpolator::sampleFromGrayscaleImage(tr, vec2(i,j));
                if(val != 0) {
                    u += val;
                    P += pow(val,2);
                    n ++;
                }
            }
        }
        u /= n;
        float sigma = sqrt((P-n*pow(u, 2))/(n-1));
        float f = desiredSigma.get()/sigma;
        for (unsigned int j = 0; j < texDims_.y - 1; j++) {
            for (unsigned int i = 0; i < texDims_.x - 1; i++) {
                float val = Interpolator::sampleFromGrayscaleImage(tr, vec2(i,j));
                if(val != 0) {
                    float new_P = desiredU.get() + f*(val - u);
                    if(new_P > 255){
                        new_P = 255;
                    }
                    if (new_P < 0){
                        new_P = 0;
                    }
                    lr->setFromDVec4(size2_t(i, j), dvec4(new_P, new_P, new_P, 255));
                }
            }
        }
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
	float pixelLength = (float) (dims.x - 1) / (texDims_.x);
	float pixelHeight = (float) (dims.y -1) / (texDims_.y);
	float stepSize = pixelLength < pixelHeight ? pixelLength : pixelHeight;
    bool ** exploredPixels = new bool*[texDims_.x];
    for (unsigned int i = 0; i < texDims_.x; i++) {
        exploredPixels[i] = new bool[texDims_.y];
        for (unsigned int j = 0; j < texDims_.y; j++)
            exploredPixels[i][j] = false;
    }
    
    for (unsigned int j = 0; j < texDims_.y; j++) {
        for (unsigned int i = 0; i < texDims_.x; i++) {
            if (exploredPixels[i][j] == true) continue;
            // Create streamline
			auto vertices = Integrator::createStreamLine(vec2(i * pixelLength, j * pixelHeight), vol.get(), 1000, stepSize, exploredPixels, pixelLength, pixelHeight);
            if(vertices.size() == 0) continue;
            float currentSum = 0;
            std::queue<float> valueQueue;
            unsigned int queueHeader = 0;
            // Init first values in queue
            while(queueHeader < kernelSize.get() / 2 && queueHeader < vertices.size()) {
                float val = Interpolator::sampleFromGrayscaleImage(tr, vec2(vertices[queueHeader].x / pixelLength, vertices[queueHeader].y / pixelHeight));
                currentSum += val;
                valueQueue.push(val);
                queueHeader++;
            }
            // find first point value
            float pointValue = currentSum / valueQueue.size();
            unsigned int pixelIdxX = vertices[0].x / pixelLength;
            unsigned int pixelIdxY = vertices[0].y / pixelHeight;
            // Draw first point
            lr->setFromDVec4(size2_t(pixelIdxX, pixelIdxY), dvec4(pointValue, pointValue, pointValue, 255));
			for(size_t l = 1; l < vertices.size(); l++) {
				vec2 vertex = vertices[l];
				pixelIdxX = vertex.x / pixelLength;
				pixelIdxY = vertex.y / pixelHeight;
                if(pixelIdxX > texDims_.x || pixelIdxY > texDims_.y || pixelIdxY < 0 || pixelIdxX < 0) continue;
                if(valueQueue.size() > kernelSize.get()) {
                    currentSum -= valueQueue.front();
                    valueQueue.pop();
                }
                if(queueHeader < vertices.size()) {
                    float val = Interpolator::sampleFromGrayscaleImage(tr, vec2(vertices[queueHeader].x / pixelLength, vertices[queueHeader].y / pixelHeight));
                    currentSum += val;
                    valueQueue.push(val);
                } else if(valueQueue.size() > kernelSize.get() / 2) {
                    currentSum -= valueQueue.front();
                    valueQueue.pop();
                }
                queueHeader++;
                pointValue = currentSum / valueQueue.size();
                lr->setFromDVec4(size2_t(pixelIdxX, pixelIdxY), dvec4(pointValue, pointValue, pointValue, 255));

			}
        }
    }
    if(useContrast.get()){
         contrast(lr, outImage -> getRepresentation<ImageRAM>());
    }
    licOut_.setData(outImage);
}

}  // namespace inviwo
