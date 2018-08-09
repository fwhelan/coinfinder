
#include <iostream>
#include <fstream>
#include <vector>

#include "bugfix.h"

#include <Python.h>

void BugFix::test1() {
	std::cerr << "In test1" << std::endl;
	 /*Embedded Python*/
        std::string returnval = "";
        Py_Initialize();
        PyObject* pValue;
	double edges_size = 0;

        PyRun_SimpleString("import sys");
	std::string path = "";
        std::string python_path = "sys.path.append(\""+path+"\")";
        PyRun_SimpleString(python_path.c_str());
        std::cerr << "Python path established." << std::endl;
        PyObject* pName = PyUnicode_DecodeFSDefault("common_ancestor");
        //PyObject* pName = PyUnicode_DecodeFSDefault("both");
                std::cerr << "Python pName established." << std::endl;
        PyObject* pModule = PyImport_Import(pName);
                std::cerr << "Python pModule established." << std::endl;
        Py_DECREF(pName);
        if (pModule != NULL) {
                PyObject* pFunc = PyObject_GetAttrString(pModule, "calc");
                if (pFunc && PyCallable_Check(pFunc)) {
                        PyObject* pArgs = PyTuple_New(1+edges_size);
                        PyTuple_SetItem(pArgs, 0, PyUnicode_FromString("((MDR102,MDR18)NODE_1,(MDR15,MDR16)NODE_2,((MDR100,MDR101)NODE_4,(MDR14,(MDR105,MDR1)NODE_6)NODE_5)NODE_3);"));
                        for (size_t i=0; i < edges_size; i++) {
                        }
                                std::cerr << "Calling common_ancestor.py..." << std::endl;
                        pValue = PyObject_CallObject(pFunc, pArgs);
                                std::cerr << "Returning from common_ancestor.py..." << std::endl;
                        Py_DECREF(pArgs);
                        if (pValue == NULL) {
                                std::cerr << "pValue is null" << std::endl;
                                PyErr_Print();
                        } else {
                                if (PyUnicode_Check(pValue) == 1) {
                                        returnval = PyUnicode_AsUTF8(pValue);
                                } else { //if (PyUnicode_Check(pValue) == 1) {
                                        PyErr_Print();
                                        throw std::logic_error( "There was an error in Python's common_ancestor. This is unexpected behaviour; please submit a bug report." );
                                }
                        }
                }
        }
        Py_Finalize();
}
