/*********************************************************************************
 *
 * Inviwo - Interactive Visualization Workshop
 * Version 0.6b
 *
 * Copyright (c) 2013-2014 Inviwo Foundation
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
 * Main file author: Rickard Englund
 *
 *********************************************************************************/

#include <modules/python3/pythonincluder.h>
#include <modules/python3/pythoninterface/pymodule.h>

#include <modules/python3/pyinviwo.h>
#include <inviwo/core/util/logcentral.h>
#include <modules/python3/pythoninterface/pyvalueparser.h>

namespace inviwo {

std::map<PyObject*,PyModule*> PyModule::instances_;

PyModule::PyModule(std::string moduleName):moduleName_(moduleName) {

    addMethod(new PyInfoMethod());
}

PyModule::~PyModule() {
    while (!methods_.empty()) {delete methods_.back(); methods_.pop_back();}
}

PyMethodDef* PyModule::getPyMethodDefs(){
    size_t N = methods_.size();
    PyMethodDef* embMethods = new PyMethodDef[N + 1];
    
    for (size_t i = 0; i < N; i++){
        embMethods[i] = { methods_[i]->getName2(), methods_[i]->getFunc(), methods_[i]->getFlags(), methods_[i]->getDesc2() };
    }
    embMethods[N] = { NULL, NULL, 0, NULL };
    return embMethods;

}

void PyModule::addMethod(PyMethod* method) {
    methods_.push_back(method);
}

const char* PyModule::getModuleName() {return moduleName_.c_str();}

void PyModule::printInfo() {
    for (unsigned int i = 0; i<methods_.size(); i++) {
        std::string msg = "print(\"";
        msg += methods_[i]->getName() + ":\t";
        msg += methods_[i]->getDesc();
        msg += "\")\n";
        PyRun_SimpleString(msg.c_str());
    }
}

std::vector<PyMethod*> PyModule::getPyMethods() {
    return methods_;
}

PyModule* PyModule::getModuleByPyObject(PyObject* obj) {
    return instances_[obj];
}

void PyModule::setPyObject(PyObject* obj)
{
    instances_[obj] = this;
}




}//namespace