#include <inviwo/core/datastructures/image/image.h>
#include <inviwo/core/datastructures/image/imagerepresentation.h>
#include <inviwo/core/datastructures/image/imagedisk.h>
#include <inviwo/core/datastructures/image/imageram.h>

namespace inviwo {

Image::Image(uvec2 dimensions, DataFormatBase format, ImageType comb) : Data2D(dimensions, format), imageType_(comb) {}

Data* Image::clone() const {
    Image* newImage = new Image(getDimension(), getDataFormat(), getImageType());
    
    //Do not copy all representations.
    //copyRepresentations(newImage);

    newImage->resize(getDimension());

    return newImage;
}

Image::~Image() {
    for (size_t i=0; i<representations_.size(); i++) {
        delete representations_[i];
    }
    representations_.clear();
}

void Image::resize(uvec2 dimensions) {
    setDimension(dimensions);
    for (size_t i=0; i<representations_.size(); i++) {
        ImageRepresentation* imageRepresentation = dynamic_cast<ImageRepresentation*>(representations_[i]) ;
        if ( imageRepresentation ) {
            imageRepresentation->resize(dimensions);
        }
    }
    setAllRepresentationsAsInvalid();
}

void Image::resizeImageRepresentations(Image* targetImage, uvec2 targetDim) {
    //TODO: check if getClassName() is necessary.
    //TODO: And also need to be tested on multiple representations_ such as ImageRAM, ImageDisk etc.,
    //TODO: optimize the code
    targetImage->resize(targetDim);
    ImageRepresentation* imageRepresentation = 0;
    ImageRepresentation* targetRepresentation = 0;
    std::vector<DataRepresentation*> &targetRepresentations = targetImage->representations_;

    if (targetRepresentations.size()) {
        for (size_t i=0; i<representations_.size(); i++) {
            imageRepresentation = dynamic_cast<ImageRepresentation*>(representations_[i]) ;        
            ivwAssert(imageRepresentation!=0, "Only image representations should be used here.");
            if (isRepresentationValid(i)) {
                size_t numberOfTargets = targetRepresentations.size();
                for (size_t j=0; j<numberOfTargets; j++) {
                    targetRepresentation = dynamic_cast<ImageRepresentation*>(targetRepresentations[j]) ;
                    ivwAssert(targetRepresentation!=0, "Only image representations should be used here.");
                    if (imageRepresentation->getClassName()==targetRepresentation->getClassName()) {
                        if (imageRepresentation->copyAndResizeImage(targetRepresentation)) {
                            targetImage->setRepresentationAsValid(j);
                            targetImage->lastValidRepresentation_ = targetRepresentations[j];
                        }
                    }
                }
            }
        }
    }
    else {
        ivwAssert(lastValidRepresentation_!=0, "Last valid representation is expected.");
        targetImage->setAllRepresentationsAsInvalid();
        targetImage->createDefaultRepresentation();
        ImageRepresentation* lastValidRepresentation = dynamic_cast<ImageRepresentation*>(lastValidRepresentation_);
        ImageRepresentation* cloneOfLastValidRepresentation = dynamic_cast<ImageRepresentation*>(lastValidRepresentation->clone());        
        targetImage->addRepresentation(cloneOfLastValidRepresentation);        

       targetImage->resize(targetDim);
       if (lastValidRepresentation->copyAndResizeImage(cloneOfLastValidRepresentation)) {
            targetImage->setRepresentationAsValid(targetImage->representations_.size()-1);
            targetImage->lastValidRepresentation_ = cloneOfLastValidRepresentation;
       }
    }
}

void Image::createDefaultRepresentation() const{
    representations_.push_back(createImageRAM(getDimension(), getImageType(), getDataFormat()));
}

} // namespace
