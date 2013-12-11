/**********************************************************************
 * Copyright (C) 2012-2013 Scientific Visualization Group - Link�ping University
 * All Rights Reserved.
 * 
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 * No part of this software may be reproduced or transmitted in any
 * form or by any means including photocopying or recording without
 * written permission of the copyright owner.
 *
 * Primary author : Alexander Johansson
 *
 **********************************************************************/

#include <inviwo/qt/widgets/properties/propertywidgetfactoryqt.h>

#include <inviwo/core/properties/boolproperty.h>
#include <inviwo/core/properties/buttonproperty.h>
#include <inviwo/core/properties/cameraproperty.h>
#include <inviwo/core/properties/compositeproperty.h>
#include <inviwo/core/properties/eventproperty.h>
#include <inviwo/core/properties/fileproperty.h>
#include <inviwo/core/properties/imageeditorproperty.h>
#include <inviwo/core/properties/matrixproperties.h>
#include <inviwo/core/properties/optionproperties.h>
#include <inviwo/core/properties/scalarproperties.h>
#include <inviwo/core/properties/stringproperty.h>
#include <inviwo/core/properties/texteditorproperty.h>
#include <inviwo/core/properties/transferfunctionproperty.h>
#include <inviwo/core/properties/vectorproperties.h>

#include <inviwo/qt/widgets/properties/boolpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/buttonpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/camerapropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/collapsivegroupboxwidgetqt.h>
#include <inviwo/qt/widgets/properties/colorpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/compositepropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/directorypropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/eventpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/filepropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatminmaxpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatmat2propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatmat3propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatmat4propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatvec2propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatvec3propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/floatvec4propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/intminmaxpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/intpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/intvec2propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/intvec3propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/intvec4propertywidgetqt.h>
#include <inviwo/qt/widgets/properties/optionpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/stringpropertywidgetqt.h>
#include <inviwo/qt/widgets/properties/texteditorwidgetqt.h>
#include <inviwo/qt/widgets/properties/imageeditorwidgetqt.h>
#include <inviwo/qt/widgets/properties/transferfunctionpropertywidgetqt.h>

#include <inviwo/qt/widgets/properties/syntaxhighlighter.h>


namespace inviwo {

PropertyWidgetFactoryQt::PropertyWidgetFactoryQt() {}

PropertyWidgetQt* PropertyWidgetFactoryQt::create(Property* property) {
    
    if (property->getSemantics()!=PropertySemantics::Default){
        if (dynamic_cast<FloatVec4Property*>(property)&& property->getSemantics() == PropertySemantics::Color) {
            return new ColorPropertyWidgetQt(static_cast<FloatVec4Property*>(property));
        }
        if (dynamic_cast<FileProperty*>(property)&& property->getSemantics() == PropertySemantics::Editor) {
             if (dynamic_cast<ImageEditorProperty*>(property))
                 return new ImageEditorWidgetQt(static_cast<FileProperty*>(property));
             else
                return new TextEditorWidgetQt(static_cast<FileProperty*>(property));
        }
		if ((dynamic_cast<FileProperty*>(property)||dynamic_cast<StringProperty*>(property))&& property->getSemantics() == PropertySemantics::Shader) {
			TextEditorWidgetQt* editor = new TextEditorWidgetQt(property);
			editor->getSyntaxHighligther()->setSyntax(GLSL);
			return editor;
		}
        if (dynamic_cast<IntVec4Property*>(property)&& property->getSemantics() == PropertySemantics::Color) {
            return new ColorPropertyWidgetQt(static_cast<IntVec4Property*>(property));
        }
        if (dynamic_cast<StringProperty*>(property)&& property->getSemantics() == PropertySemantics::Editor) {
            return new TextEditorWidgetQt(static_cast<StringProperty*>(property));
        }
    }
    if (dynamic_cast<BoolProperty*>(property))
        return new BoolPropertyWidgetQt(static_cast<BoolProperty*>(property));
    if (dynamic_cast<ButtonProperty*>(property))
        return new ButtonPropertyWidgetQt(static_cast<ButtonProperty*>(property));
    if (dynamic_cast<CompositeProperty*>(property))
        return new CompositePropertyWidgetQt(static_cast<CompositeProperty*>(property));
    if (dynamic_cast<DirectoryProperty*>(property))
        return new DirectoryPropertyWidgetQt(static_cast<DirectoryProperty*>(property));
    if (dynamic_cast<EventProperty*>(property))
        return new EventPropertyWidgetQt(static_cast<EventProperty*>(property));
    if (dynamic_cast<FileProperty*>(property))
        return new FilePropertyWidgetQt(static_cast<FileProperty*>(property));
	if (dynamic_cast<FloatMat2Property*>(property))
		return new FloatMat2PropertyWidgetQt(static_cast<FloatMat2Property*>(property));
    if (dynamic_cast<FloatMat3Property*>(property))
        return new FloatMat3PropertyWidgetQt(static_cast<FloatMat3Property*>(property));
    if (dynamic_cast<FloatMat4Property*>(property))
        return new FloatMat4PropertyWidgetQt(static_cast<FloatMat4Property*>(property));
    if (dynamic_cast<FloatProperty*>(property))
        return new FloatPropertyWidgetQt(static_cast<FloatProperty*>(property));
    if (dynamic_cast<FloatMinMaxProperty*>(property))
        return new FloatMinMaxPropertyWidgetQt(static_cast<FloatMinMaxProperty*>(property));
    if (dynamic_cast<FloatVec2Property*>(property))
        return new FloatVec2PropertyWidgetQt(static_cast<FloatVec2Property*>(property));
    if (dynamic_cast<FloatVec3Property*>(property))
        return new FloatVec3PropertyWidgetQt(static_cast<FloatVec3Property*>(property));
    if (dynamic_cast<FloatVec4Property*>(property))
        return new FloatVec4PropertyWidgetQt(static_cast<FloatVec4Property*>(property));
    if (dynamic_cast<IntMinMaxProperty*>(property))
        return new IntMinMaxPropertyWidgetQt(static_cast<IntMinMaxProperty*>(property));
    if (dynamic_cast<IntVec2Property*>(property))
        return new IntVec2PropertyWidgetQt(static_cast<IntVec2Property*>(property));
    if (dynamic_cast<IntVec3Property*>(property))
        return new IntVec3PropertyWidgetQt(static_cast<IntVec3Property*>(property));
    if (dynamic_cast<IntVec4Property*>(property))
        return new IntVec4PropertyWidgetQt(static_cast<IntVec4Property*>(property));
    if (dynamic_cast<IntProperty*>(property))
        return new IntPropertyWidgetQt(static_cast<IntProperty*>(property));
    if (dynamic_cast<BaseOptionProperty*>(property))
        return new OptionPropertyWidgetQt(static_cast<BaseOptionProperty*>(property));
    if (dynamic_cast<StringProperty*>(property))
        return new StringPropertyWidgetQt(static_cast<StringProperty*>(property));
    if (dynamic_cast<TransferFunctionProperty*>(property))
        return new TransferFunctionPropertyWidgetQt(static_cast<TransferFunctionProperty*>(property));

    LogWarn("No widget for property " + property->getIdentifier() + " found.")
        return 0;
}

} // namespace
