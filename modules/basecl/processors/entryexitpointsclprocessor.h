/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2014 Inviwo Foundation
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
 * Contact: Daniel J�nsson
 *
 *********************************************************************************/

#ifndef IVW_ENTRYEXITPOINTSCLPROCESSOR_H
#define IVW_ENTRYEXITPOINTSCLPROCESSOR_H

#include <modules/basecl/baseclmoduledefine.h>
#include <inviwo/core/common/inviwo.h>
#include <inviwo/core/interaction/cameratrackball.h>
#include <inviwo/core/properties/cameraproperty.h>
#include <inviwo/core/ports/geometryport.h>
#include <inviwo/core/ports/imageport.h>
#include <inviwo/core/ports/volumeport.h>
#include <inviwo/core/common/inviwoapplication.h>
#include <modules/opencl/inviwoopencl.h>
#include <modules/opencl/kernelowner.h>
#include <modules/basecl/meshentryexitpointscl.h>

namespace inviwo {

class IVW_MODULE_BASECL_API EntryExitPointsCLProcessor : public Processor, public ProcessorKernelOwner {
public:
    EntryExitPointsCLProcessor();
    virtual ~EntryExitPointsCLProcessor();

    InviwoProcessorInfo();

    virtual void initialize();
    virtual void deinitialize();

protected:
    virtual void process();



private:
    GeometryInport geometryPort_;
    ImageOutport entryPort_;
    ImageOutport exitPort_;

    CameraProperty camera_;
    IntVec2Property workGroupSize_;
    BoolProperty useGLSharing_;

    CameraTrackball trackball_;

    MeshEntryExitPointsCL entryExitPoints_;
};

}  // namespace

#endif  // IVW_ENTRYEXITPOINTSCLPROCESSOR_H