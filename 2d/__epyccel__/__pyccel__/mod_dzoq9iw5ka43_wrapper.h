#ifndef MOD_DZOQ9IW5KA43_WRAPPER_H
#define MOD_DZOQ9IW5KA43_WRAPPER_H

#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include "cwrapper_ndarrays.h"


#ifdef MOD_DZOQ9IW5KA43_WRAPPER

void bind_c_assemble_stiffness_1d(int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t);

#else

static void** Pymod_dzoq9iw5ka43_API;


/*........................................*/
static int mod_dzoq9iw5ka43_import(void)
{
    PyObject* current_path;
    PyObject* stash_path;
    current_path = PySys_GetObject("path");
    stash_path = PyList_GetItem(current_path, 0);
    Py_INCREF(stash_path);
    PyList_SetItem(current_path, 0, PyUnicode_FromString("/home/rifqui/Desktop/Public/IGASIMPLINES/2d/__epyccel__"));
    Pymod_dzoq9iw5ka43_API = (void**)PyCapsule_Import("mod_dzoq9iw5ka43._C_API", 0);
    PyList_SetItem(current_path, 0, stash_path);
    return Pymod_dzoq9iw5ka43_API != NULL ? 0 : -1;
}
/*........................................*/

#endif
#endif // MOD_DZOQ9IW5KA43_WRAPPER_H
