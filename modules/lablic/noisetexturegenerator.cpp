/*********************************************************************
 *  Author  : Himangshu Saikia
 *  Init    : Monday, October 02, 2017 - 13:31:36
 *
 *  Project : KTH Inviwo Modules
 *
 *  License : Follows the Inviwo BSD license model
 *********************************************************************
 */

#include <lablic/noisetexturegenerator.h>
#include <inviwo/core/datastructures/image/layerram.h>
#include <inviwo/core/util/utilities.h>

namespace inviwo {

// The Class Identifier has to be globally unique. Use a reverse DNS naming scheme
const ProcessorInfo NoiseTextureGenerator::processorInfo_{
    "org.inviwo.NoiseTextureGenerator",  // Class identifier
    "Noise Texture Generator",           // Display name
    "KTH Labs",                                  // Category
    CodeState::Experimental,             // Code state
    Tags::None,                          // Tags
};

const ProcessorInfo NoiseTextureGenerator::getProcessorInfo() const { return processorInfo_; }

NoiseTextureGenerator::NoiseTextureGenerator()
	: Processor()
	, texOut_("texOut")
	, texSize_("texSize", "Texture Size", vec2(512, 512), vec2(1, 1), vec2(2048, 2048), vec2(1, 1))
	, noise_mode("noise_mode", "Noise Mode")
	, seed_random("seed_random", "Use random seed", false)
	, random_seed("random_seed", "Random Seed", 500, 0, 1000, 1)

    // TODO: Register additional properties

{
    // Register ports
    addPort(texOut_);

    // Register properties
	noise_mode.addOption("grayscale", "Grayscale", 0);
	noise_mode.addOption("black-white", "Black-While", 1);
	addProperty(noise_mode);
    addProperty(texSize_);
	addProperty(seed_random);
	addProperty(random_seed);
	
	//Only display relevant properties
	util::hide(random_seed);
	seed_random.onChange([this]() {
		if (seed_random.get() == false) {
			util::hide(random_seed);
		}
		else {
			util::show(random_seed);
		}
	});
}

void NoiseTextureGenerator::process() {
    // The output of the generation process is an Image
    auto outImage = std::make_shared<Image>();

    // In Inviwo, images consist of multiple layers, but only the
    // color layer is relevant for us
    auto outLayer = outImage->getColorLayer();

    outLayer->setDimensions(size2_t(texSize_.get().x, texSize_.get().y));
    // With the data format DataVec4UInt8 values for RGB-alpha range between 0 and 255
    outLayer->setDataFormat(DataVec4UInt8::get());

    // Just like we did with the volume in other assignments we need to retrieve an editable
    // representation of the object we want to modify (here a layer)
    auto lr = outLayer->getEditableRepresentation<LayerRAM>();

	//Seed random number generator
	if (seed_random.get() == true)
		srand(random_seed.get());
	else
		srand(time(NULL));

	//Color the pixels
    for (int j = 0; j < texSize_.get().y; j++) {
        for (int i = 0; i < texSize_.get().x; i++) {

            // Randomly sample values for the texture			
			float r_temp = ((float)rand()) / (float)RAND_MAX;
			float random = r_temp * 255;

			float value = random;
			if (noise_mode.get() == 1)
				value = r_temp > 0.5 ? 0.0f : 255.0f;

			// A value within the ouput image is set by specifying pixel position and color
			lr->setFromDVec4(size2_t(i, j), dvec4(value, value, value, 255));
        }
    }

    texOut_.setData(outImage);
}

}  // namespace inviwo
