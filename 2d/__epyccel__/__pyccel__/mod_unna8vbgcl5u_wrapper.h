#ifndef MOD_UNNA8VBGCL5U_WRAPPER_H
#define MOD_UNNA8VBGCL5U_WRAPPER_H

#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include "cwrapper_ndarrays.h"


#ifdef MOD_UNNA8VBGCL5U_WRAPPER

void bind_c_assemble_vector_ex001(int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t);

#else

static void** Pymod_unna8vbgcl5u_API;


/*........................................*/
static int mod_unna8vbgcl5u_import(void)
{
    PyObject* current_path;
    PyObject* stash_path;
    current_path = PySys_GetObject("path");
    stash_path = PyList_GetItem(current_path, 0);
    Py_INCREF(stash_path);
    PyList_SetItem(current_path, 0, PyUnicode_FromString("/home/rifqui/Desktop/Public/IGASIMPLINES/2d/__epyccel__"));
    Pymod_unna8vbgcl5u_API = (void**)PyCapsule_Import("mod_unna8vbgcl5u._C_API", 0);
    PyList_SetItem(current_path, 0, stash_path);
    return Pymod_unna8vbgcl5u_API != NULL ? 0 : -1;
}
/*........................................*/

#endif
#endif // MOD_UNNA8VBGCL5U_WRAPPER_H
