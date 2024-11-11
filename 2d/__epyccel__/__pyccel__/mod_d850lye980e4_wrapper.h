#ifndef MOD_D850LYE980E4_WRAPPER_H
#define MOD_D850LYE980E4_WRAPPER_H

#include "numpy_version.h"
#include "numpy/arrayobject.h"
#include "cwrapper.h"
#include "cwrapper_ndarrays.h"


#ifdef MOD_D850LYE980E4_WRAPPER

void bind_c_assemble_vector_ex001(int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t, void*, int64_t, int64_t, int64_t, int64_t);

#else

static void** Pymod_d850lye980e4_API;


/*........................................*/
static int mod_d850lye980e4_import(void)
{
    PyObject* current_path;
    PyObject* stash_path;
    current_path = PySys_GetObject("path");
    stash_path = PyList_GetItem(current_path, 0);
    Py_INCREF(stash_path);
    PyList_SetItem(current_path, 0, PyUnicode_FromString("/home/rifqui/Desktop/Public/IGASIMPLINES/2d/__epyccel__"));
    Pymod_d850lye980e4_API = (void**)PyCapsule_Import("mod_d850lye980e4._C_API", 0);
    PyList_SetItem(current_path, 0, stash_path);
    return Pymod_d850lye980e4_API != NULL ? 0 : -1;
}
/*........................................*/

#endif
#endif // MOD_D850LYE980E4_WRAPPER_H
