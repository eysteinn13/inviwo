/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 *
 * Copyright (c) 2015-2016 Inviwo Foundation
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *********************************************************************************/

#ifndef IVW_VOLUMEBASISPROPERTY_H
#define IVW_VOLUMEBASISPROPERTY_H

#include <modules/base/basemoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/properties/compositeproperty.h>
#include <inviwo/core/properties/ordinalproperty.h>
#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/datastructures/spatialdata.h>

namespace inviwo {

/**
 * \ingroup properties
 * A CompositeProperty holding the properties needed to represent a bases matrix.
 */
class IVW_MODULE_BASE_API BasisProperty : public CompositeProperty {
public:
    InviwoPropertyInfo();
    BasisProperty(std::string identifier, std::string displayName,
                        InvalidationLevel invalidationLevel = InvalidationLevel::InvalidResources,
                        PropertySemantics semantics = PropertySemantics::Default);
    BasisProperty(const BasisProperty& rhs);
    BasisProperty& operator=(const BasisProperty& that);
    virtual BasisProperty* clone() const override;
    virtual ~BasisProperty() = default;

    void updateForNewEntity(const SpatialEntity<3>& volume, bool deserialize = false);

    void updateEntity(SpatialEntity<3>& volume);

    BoolProperty overRideDefaults_;
    FloatVec3Property a_;
    FloatVec3Property b_;
    FloatVec3Property c_;
    FloatVec3Property offset_;

    FloatVec3Property overrideA_;
    FloatVec3Property overrideB_;
    FloatVec3Property overrideC_;
    FloatVec3Property overrideOffset_;

    mat4 getBasisAndOffset() const;

private:
    void onOverrideChange();
};

} // namespace

#endif // IVW_VOLUMEBASISPROPERTY_H

