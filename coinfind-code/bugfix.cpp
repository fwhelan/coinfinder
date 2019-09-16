
#include <iostream>
#include <fstream>
#include <vector>

#include "bugfix.h"

#include <Python.h>

void BugFix::test1() {
	std::cerr << "In test1, calling first python piece" << std::endl;
	 /*Embedded Python*/
        Py_Initialize();
        PyObject* pValue;
        
        PyRun_SimpleString("import sys");
	std::string path = "";
        std::string python_path = "sys.path.append(\""+path+"\")";
        PyRun_SimpleString(python_path.c_str());
        
        PyObject* pName = PyUnicode_DecodeFSDefault("phylomax");
        //PyObject* pName = PyUnicode_DecodeFSDefault("both");
        PyObject* pModule = PyImport_Import(pName);
        Py_DECREF(pName);
        if (pModule != NULL) {
                PyObject* pFunc = PyObject_GetAttrString(pModule, "calc");
                if (pFunc && PyCallable_Check(pFunc)) {
                        PyObject* pArgs = PyTuple_New(1);
                        PyTuple_SetItem(pArgs, 0, PyUnicode_FromString("((MDR102,MDR18)NODE_1,(MDR15,MDR16)NODE_2,((MDR100,MDR101)NODE_4,(MDR14,(MDR105,MDR1)NODE_6)NODE_5)NODE_3);"));
                        pValue = PyObject_CallObject(pFunc, pArgs);
                        Py_DECREF(pArgs);
                        if (pValue == NULL) {
                                std::cerr << "pValue is null" << std::endl;
                                PyErr_Print();
                        } else {
                                if (PyUnicode_Check(pValue) == 1) { //return value is a string; there was an error
                                        PyErr_Print();
                                        throw std::logic_error( "Error: python's phylomax returned a Unicode. This is unexpected behaviour; please submit a bug report." );
                                } else if (PyFloat_Check(pValue) == 1) {
                                        PyErr_Print();
                                        throw std::logic_error( "Error: python's phylomax returned a Float. This is unexpected behaviour; please submit a bug report." );
                                } else if (PyList_Check(pValue) == 1) { //return value is a list
                                        /*Save returned list into a map with beta names as paired keys and dist as value*/
                                        for(int a=0; a<PyList_Size(pValue); a=a+3) {
                                                PyObject *value1 = PyList_GetItem(pValue, a);
                                                PyObject *value2 = PyList_GetItem(pValue, a+1);
                                                PyObject *value3 = PyList_GetItem(pValue, a+2);
                                                //phylogenetic_distances[PyFloat_AsDouble(value3)] = std::make_pair(PyUnicode_AsUTF8(value1),PyUnicode_AsUTF8(value2));
                                        }
                                
                                } else {
                                        throw std::logic_error( "Error: a list isn't being returned from phylomax. This is unexpected behaviour; please submit a bug report." );
                                }
                                Py_DECREF(pValue);
                        }
                }
        } else {
                PyErr_Print();
                throw std::logic_error( "There was an error in Python's phylomax. This is unexpected behaviour; please submit a bug report." );
        }
        Py_Finalize();


	std::cerr << "Calling second python piece" << std::endl;

	 /*Embedded Python*/
        std::string returnval = "";
        Py_Initialize();
        // PyObject* pValue;
	double edges_size = 0;

        PyRun_SimpleString("import sys");
        PyRun_SimpleString(python_path.c_str());
        std::cerr << "Python path established." << std::endl;
        pName = PyUnicode_DecodeFSDefault("common_ancestor");
        //PyObject* pName = PyUnicode_DecodeFSDefault("both");
                std::cerr << "Python pName established." << std::endl;
        pModule = PyImport_Import(pName);
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
