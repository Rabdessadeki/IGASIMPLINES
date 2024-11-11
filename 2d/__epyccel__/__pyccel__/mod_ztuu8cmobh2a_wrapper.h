#ifndef MOD_ZTUU8CMOBH2A_WRAPPER_H
#define MOD_ZTUU8CMOBH2A_WRAPPER_H

#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include "cwrapper_ndarrays.h"


#ifdef MOD_ZTUU8CMOBH2A_WRAPPER

void bind_c_assemble_stiffness_1d(int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t);

#else

static void** Pymod_ztuu8cmobh2a_API;


/*........................................*/
static int mod_ztuu8cmobh2a_import(void)
{
    PyObject* current_path;
    PyObject* stash_path;
    current_path = PySys_GetObject("path");
    stash_path = PyList_GetItem(current_path, 0);
    Py_INCREF(stash_path);
    PyList_SetItem(current_path, 0, PyUnicode_FromString("/home/rifqui/Desktop/Public/IGASIMPLINES/2d/__epyccel__"));
    Pymod_ztuu8cmobh2a_API = (void**)PyCapsule_Import("mod_ztuu8cmobh2a._C_API", 0);
    PyList_SetItem(current_path, 0, stash_path);
    return Pymod_ztuu8cmobh2a_API != NULL ? 0 : -1;
}
/*........................................*/

#endif
#endif // MOD_ZTUU8CMOBH2A_WRAPPER_H
