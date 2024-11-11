#define PY_ARRAY_UNIQUE_SYMBOL CWRAPPER_ARRAY_API
#define MOD_4J4ZIFXXOQJF_WRAPPER

#include "mod_4j4zifxxoqjf_wrapper.h"
#include <stdlib.h>
#include <stdint.h>
#include "ndarrays.h"


/*........................................*/


/*........................................*/

/*........................................*/
static PyObject* bind_c_assemble_stiffness_1d_wrapper(PyObject* self, PyObject* args, PyObject* kwargs)
{
    PyObject* ne_obj;
    PyObject* degree_obj;
    PyObject* spans_obj;
    PyObject* basis_obj;
    PyObject* weights_obj;
    PyObject* points_obj;
    PyObject* matrix_obj;
    int64_t ne;
    int64_t degree;
    t_ndarray spans = {.shape = NULL};
    void* bound_spans;
    int64_t bound_spans_shape_1;
    int64_t bound_spans_stride_1;
    t_ndarray basis = {.shape = NULL};
    void* bound_basis;
    int64_t bound_basis_shape_1;
    int64_t bound_basis_shape_2;
    int64_t bound_basis_shape_3;
    int64_t bound_basis_shape_4;
    int64_t bound_basis_stride_1;
    int64_t bound_basis_stride_2;
    int64_t bound_basis_stride_3;
    int64_t bound_basis_stride_4;
    t_ndarray weights = {.shape = NULL};
    void* bound_weights;
    int64_t bound_weights_shape_1;
    int64_t bound_weights_shape_2;
    int64_t bound_weights_stride_1;
    int64_t bound_weights_stride_2;
    t_ndarray points = {.shape = NULL};
    void* bound_points;
    int64_t bound_points_shape_1;
    int64_t bound_points_shape_2;
    int64_t bound_points_stride_1;
    int64_t bound_points_stride_2;
    t_ndarray matrix = {.shape = NULL};
    void* bound_matrix;
    int64_t bound_matrix_shape_1;
    int64_t bound_matrix_shape_2;
    int64_t bound_matrix_stride_1;
    int64_t bound_matrix_stride_2;
    static char *kwlist[] = {
        "ne",
        "degree",
        "spans",
        "basis",
        "weights",
        "points",
        "matrix",
        NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOO", kwlist, &ne_obj, &degree_obj, &spans_obj, &basis_obj, &weights_obj, &points_obj, &matrix_obj))
    {
        return NULL;
    }
    if (PyIs_NativeInt(ne_obj))
    {
        ne = PyInt64_to_Int64(ne_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument ne");
        return NULL;
    }
    if (PyIs_NativeInt(degree_obj))
    {
        degree = PyInt64_to_Int64(degree_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument degree");
        return NULL;
    }
    if (pyarray_check(spans_obj, NPY_LONG, INT64_C(1), NO_ORDER_CHECK))
    {
        spans = pyarray_to_ndarray(spans_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type int for argument spans");
        return NULL;
    }
    bound_spans = nd_data(&spans);
    bound_spans_shape_1 = nd_ndim(&spans, INT64_C(0));
    bound_spans_stride_1 = nd_nstep_F(&spans, INT64_C(0));
    if (pyarray_check(basis_obj, NPY_DOUBLE, INT64_C(4), NPY_ARRAY_C_CONTIGUOUS))
    {
        basis = pyarray_to_ndarray(basis_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument basis");
        return NULL;
    }
    bound_basis = nd_data(&basis);
    bound_basis_shape_1 = nd_ndim(&basis, INT64_C(0));
    bound_basis_shape_2 = nd_ndim(&basis, INT64_C(1));
    bound_basis_shape_3 = nd_ndim(&basis, INT64_C(2));
    bound_basis_shape_4 = nd_ndim(&basis, INT64_C(3));
    bound_basis_stride_1 = nd_nstep_C(&basis, INT64_C(0));
    bound_basis_stride_2 = nd_nstep_C(&basis, INT64_C(1));
    bound_basis_stride_3 = nd_nstep_C(&basis, INT64_C(2));
    bound_basis_stride_4 = nd_nstep_C(&basis, INT64_C(3));
    if (pyarray_check(weights_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        weights = pyarray_to_ndarray(weights_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument weights");
        return NULL;
    }
    bound_weights = nd_data(&weights);
    bound_weights_shape_1 = nd_ndim(&weights, INT64_C(0));
    bound_weights_shape_2 = nd_ndim(&weights, INT64_C(1));
    bound_weights_stride_1 = nd_nstep_C(&weights, INT64_C(0));
    bound_weights_stride_2 = nd_nstep_C(&weights, INT64_C(1));
    if (pyarray_check(points_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        points = pyarray_to_ndarray(points_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument points");
        return NULL;
    }
    bound_points = nd_data(&points);
    bound_points_shape_1 = nd_ndim(&points, INT64_C(0));
    bound_points_shape_2 = nd_ndim(&points, INT64_C(1));
    bound_points_stride_1 = nd_nstep_C(&points, INT64_C(0));
    bound_points_stride_2 = nd_nstep_C(&points, INT64_C(1));
    if (pyarray_check(matrix_obj, NPY_DOUBLE, INT64_C(2), NPY_ARRAY_C_CONTIGUOUS))
    {
        matrix = pyarray_to_ndarray(matrix_obj);
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "Expected an argument of type float for argument matrix");
        return NULL;
    }
    bound_matrix = nd_data(&matrix);
    bound_matrix_shape_1 = nd_ndim(&matrix, INT64_C(0));
    bound_matrix_shape_2 = nd_ndim(&matrix, INT64_C(1));
    bound_matrix_stride_1 = nd_nstep_C(&matrix, INT64_C(0));
    bound_matrix_stride_2 = nd_nstep_C(&matrix, INT64_C(1));
    bind_c_assemble_stiffness_1d(ne, degree, bound_spans, bound_spans_shape_1, bound_spans_stride_1, bound_basis, bound_basis_shape_1, bound_basis_shape_2, bound_basis_shape_3, bound_basis_shape_4, bound_basis_stride_1, bound_basis_stride_2, bound_basis_stride_3, bound_basis_stride_4, bound_weights, bound_weights_shape_1, bound_weights_shape_2, bound_weights_stride_1, bound_weights_stride_2, bound_points, bound_points_shape_1, bound_points_shape_2, bound_points_stride_1, bound_points_stride_2, bound_matrix, bound_matrix_shape_1, bound_matrix_shape_2, bound_matrix_stride_1, bound_matrix_stride_2);
    free_pointer(&spans);
    free_pointer(&basis);
    free_pointer(&weights);
    free_pointer(&points);
    free_pointer(&matrix);
    Py_INCREF(Py_None);
    return Py_None;
}
/*........................................*/

/*........................................*/

static PyMethodDef mod_4j4zifxxoqjf_methods[] = {
    {
        "assemble_stiffness_1D",
        (PyCFunction)bind_c_assemble_stiffness_1d_wrapper,
        METH_VARARGS | METH_KEYWORDS,
        ""
    },
    { NULL, NULL, 0, NULL}
};

/*........................................*/

static struct PyModuleDef mod_4j4zifxxoqjf_module = {
    PyModuleDef_HEAD_INIT,
    /* name of module */
    "mod_4j4zifxxoqjf",
    /* module documentation, may be NULL */
    NULL,
    /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    0,
    mod_4j4zifxxoqjf_methods,
};

/*........................................*/

PyMODINIT_FUNC PyInit_mod_4j4zifxxoqjf(void)
{
    PyObject* mod;
    static void* Pymod_4j4zifxxoqjf_API[0];
    PyObject* c_api_object_0001;
    mod = PyModule_Create(&mod_4j4zifxxoqjf_module);
    if (mod == NULL)
    {
        return NULL;
    }
    c_api_object_0001 = PyCapsule_New((void *)Pymod_4j4zifxxoqjf_API, "mod_4j4zifxxoqjf._C_API", NULL);
    if (PyModule_AddObject(mod, "_C_API", c_api_object_0001) < INT64_C(0))
    {
        Py_DECREF(mod);
        return NULL;
    }
    Py_INCREF(c_api_object_0001);
    import_array();
    return mod;
}
