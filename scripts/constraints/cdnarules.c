#include <Python.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION
#include <numpy/arrayobject.h>   /* NumPy  as seen from C */
//#if defined(__MACH__)
//    #include <stdlib.h>
//#else
//    #include <malloc.h>
//#endif
#define T uint64_t
#define BYTE uint8_t // unsigned char

static PyObject* bitSet(PyObject* self, PyObject *args) {
   T v;
   unsigned int b;
   if (!PyArg_ParseTuple(args, "Ki", &v, &b)) {
      return NULL;
   }
   bool c = ((v >> b) & 1) == 1;
   PyObject *return_val = Py_BuildValue("b", c);
   return return_val;
}

static inline int bitsSet_internal(T v) {
   v = v - ((v >> 1) & (T)~(T)0/3);                           // temp
   v = (v & (T)~(T)0/15*3) + ((v >> 2) & (T)~(T)0/15*3);      // temp
   v = (v + (v >> 4)) & (T)~(T)0/255*15;                      // temp
   int c = (T)(v * ((T)~(T)0/255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
   return c;
}

static PyObject* bitsSet(PyObject* self, PyObject *args) {
   T v;
   if (!PyArg_ParseTuple(args, "K", &v)) {
      return NULL;
   }
   int c = bitsSet_internal(v);
   PyObject *return_val = Py_BuildValue("i", c);
   return return_val;
}

static PyObject* microsatellite(PyObject* self,  PyObject *args)
{
   char *text;
   int lengthToLookFor;
   if (!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor)) {
      return NULL;
   }
   int i = 0;
   int n = strlen(text);
   int res = 1;
   char *resChars = PyMem_RawMalloc((lengthToLookFor + 1) * sizeof(char));
   if (resChars == NULL)
       return PyErr_NoMemory();
   strncpy( resChars, text, lengthToLookFor );
   resChars[lengthToLookFor] = '\0';
   //text[:lengthToLookFor];
   int maxLenght = 0;
   int ret;
   while (i <= n - 2 * lengthToLookFor) {
       ret = memcmp(&text[i],&text[i+lengthToLookFor], lengthToLookFor);
       if (ret == 0) {
           // we found one
            res += 1;
        }else{
            if (maxLenght < res) {
                maxLenght = res;
                //resChars = text[i: i + lengthToLookFor];
                strncpy( resChars, &text[i], lengthToLookFor );
            }
            res = 1;
        }
        i += lengthToLookFor;
   }
   if (maxLenght < res) {
        maxLenght = res;
        //resChars = text[i: i + lengthToLookFor];
        strncpy( resChars, &text[i], lengthToLookFor );
    }
   PyObject *return_val = Py_BuildValue("(is)", maxLenght, resChars);
   PyMem_RawFree(resChars);
   return return_val;
}

static PyObject* longestSequenceOfChar(PyObject* self,  PyObject *args)
{
   char *text; //dna_seq
   char *char_x; // dna base to check (* for all)
   if (!PyArg_ParseTuple(args, "ss", &text, &char_x)) {
      return NULL;
   }
   if (strlen(char_x) != 1) {
      return NULL;
   }
   int c = 0; // max homopolymer length found so far
   char res = char_x[0];  // text[0]
   int n = strlen(text);
   int curr = 1;
   int i = 0;
   while (i < n-1) {
        if (i < n - 1 && text[i] == text[i+1]) {
           curr += 1;
        } else {
            if (curr > c) {
                if (char_x[0] == '*' || text[i] == char_x[0]) {
                    c = curr;
                    res = text[i];
                }
            }
            curr = 1;
        }
        i += 1;
   }
   if (curr > c && (char_x[0] == '*' || text[i] == char_x[0])) {
       c = curr;
       res = text[i];
   }
   PyObject *return_val = Py_BuildValue("(ci)", res, c);
   return return_val;
}

static PyObject* repeatRegion(PyObject* self,  PyObject *args)
{
   char *text;
   int res = 0;
   int lengthToLookFor;
   if (!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor)) {
      return NULL;
   }
   int len = strlen(text);
   char *subseq = PyMem_RawMalloc((lengthToLookFor+1) * sizeof(char));
   if (subseq == NULL)
       return PyErr_NoMemory();
   strncpy( subseq, text, lengthToLookFor );
   subseq[lengthToLookFor] = '\0';
   for (int i = 0; i < len-lengthToLookFor;i++) {
        strncpy( subseq, &text[i], lengthToLookFor);
        if (strstr(&text[i+1], subseq) != NULL) {
            res = 1;
            break;
        }
   }
   PyObject *return_val = Py_BuildValue("i", res);
   PyMem_RawFree(subseq);
   return return_val;
}

static PyObject* smallRepeatRegion(PyObject* self,  PyObject *args)
{
   char *text;
   float res = 1.0;
   int lengthToLookFor;
   if (!PyArg_ParseTuple(args, "si", &text, &lengthToLookFor)) {
      return NULL;
   }
   int len = strlen(text);
   char *subseq = PyMem_RawMalloc((lengthToLookFor+1) * sizeof(char));
   if (subseq == NULL)
       return PyErr_NoMemory();
   strncpy(subseq, text, lengthToLookFor );
   subseq[lengthToLookFor] = '\0';
   for (int i = 0; i <= len-lengthToLookFor;i++) {
        strncpy(subseq, &text[i], lengthToLookFor);
        if (strstr(&text[i+1], subseq) != NULL) {
            res += 1.0;
        }
   }
   if (res * lengthToLookFor / strlen(text) > 0.44) {
       res = 1.0;
   } else {
       res = res * lengthToLookFor / strlen(text) * 0.5;
   }
   PyObject *return_val = Py_BuildValue("d", res);
   PyMem_RawFree(subseq);
   return return_val;
}

static PyObject* getQUAT(PyObject* self,  PyObject *args)
{
   int bit1, bit2;
   char *res = PyMem_RawMalloc(2 * sizeof(char));
   if (res == NULL)
       return PyErr_NoMemory();
   res[0] = 'A';
   res[1] = '\0';
   if (!PyArg_ParseTuple(args, "pp", &bit1, &bit2)) {
      return NULL;
   }
   if (bit1 && bit2 ) {
       res[0] = 'T';
   } else if (!bit1 && bit2) {
       res[0] = 'C';
   } else if (bit1 && !bit2) {
       res[0] = 'G';
   }
   PyObject *return_val = Py_BuildValue("s", res);
   PyMem_RawFree(res);
   return return_val;
}

static PyObject* byte2QUATS(PyObject* self,  PyObject *args)
{
   int byte;
   int bit1, bit2;
   char *res = PyMem_RawMalloc(5 * sizeof(char));
   if (res == NULL)
       return PyErr_NoMemory();
   res[0] = 'A';
   res[1] = 'A';
   res[2] = 'A';
   res[3] = 'A';
   res[4] = '\0';
   if (!PyArg_ParseTuple(args, "i", &byte)) {
      return NULL;
   }
   int pos = 6;
   for (int i=0; i < 4; i++) {
        bit1 = ((byte >> (pos+1))  & 0x01);
        bit2 = ((byte >> pos)  & 0x01);
        if (bit1 && bit2 ) {
           res[i] = 'T';
        } else if (!bit1 && bit2) {
           res[i] = 'C';
        } else if (bit1 && !bit2) {
           res[i] = 'G';
        }
        pos -= 2;
   }
   PyObject *return_val = Py_BuildValue("s", res);
   PyMem_RawFree(res);
   return return_val;
}

static PyObject* strContainsSub(PyObject* self,  PyObject *args)
{
   char *text;
   char *substr;
   int res = 0;
   if (!PyArg_ParseTuple(args, "ss", &text, &substr)) {
      return NULL;
   }
   if (strstr(text, substr) != NULL) {
      res = 1;
   }
   PyObject *return_val = Py_BuildValue("O", res ? Py_True: Py_False);
   return return_val;
}

static char cdnarules_sat_docs[] =
   "microsatellite(text, lengthToLookFor): Finds the maximum microsatellite with the given lenght!\n";
static char cdnarules_lseq_docs[] =
   "longestSequenceOfChar(text, character_to_look_for): Finds the maximum sequence of a char (or any char if given '*')!\n";
static char cdnarules_rreg_docs[] =
   "repeatRegion(text, lengthToLookFor): returns true if a region of 'lengthToLookFor' is repeated within text (includes overlapping texts)";
static char cdnarules_srreg_docs[] =
   "smallRepeatRegion(text, lengthToLookFor): returns a error-value based on the number repeats the text contains";
static char cdnarules_getq_docs[] =
   "getQUAT(bit1,bit2): returns the DNA base for the given bits";
static char cdnarules_byte2quats_docs[] =
   "byte2QUATS(byte): Converts a given byte to DNA representation";
static char cdnarules_strcsubstr[] =
   "strContainsSub(text,substr): Returns true if the substr is present in text";
static char bitsSet_docs[] =
    "bitsSet(integer): returns the number of bits set int given integer";
static char bitSet_docs[] =
"bitSet(X,b): returns if bit b is set in X";
static PyMethodDef cdnarules_funcs[] = {
   //{"microsatellite", (PyCFunction)microsatellite, METH_NOARGS, cdnarules_docs},
   {"bitsSet", bitsSet, METH_VARARGS, bitsSet_docs},
   {"microsatellite", microsatellite, METH_VARARGS, cdnarules_sat_docs},
   {"longestSequenceOfChar", longestSequenceOfChar, METH_VARARGS, cdnarules_lseq_docs},
   {"repeatRegion", repeatRegion, METH_VARARGS, cdnarules_rreg_docs},
   {"smallRepeatRegion", smallRepeatRegion, METH_VARARGS, cdnarules_srreg_docs},
   {"getQUAT", getQUAT, METH_VARARGS, cdnarules_getq_docs},
   {"byte2QUATS", byte2QUATS, METH_VARARGS, cdnarules_byte2quats_docs},
   {"strContainsSub", strContainsSub, METH_VARARGS, cdnarules_strcsubstr},
   {"bitsSet", bitsSet, METH_VARARGS, bitsSet_docs},
   {"bitSet", bitSet, METH_VARARGS, bitSet_docs},
   {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cdnarules =
{
    PyModuleDef_HEAD_INIT,
    "cdnarules", /* name of module */
    "Extension module for fast DNARules processing!",
    //"usage: Combinations.uniqueCombinations(lstSortableItems, comboSize)\n", /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    cdnarules_funcs
};

PyMODINIT_FUNC PyInit_cdnarules(void)
{
    //PyObject* module =  PyModule_Create(&cdnarules);
    import_array();
    //return module;
    return PyModule_Create(&cdnarules);
}

/*void initcdnarules(void)
{
   Py_Initialize3("cdnarules", cdnarules_funcs, "Extension module for fast DNARules processing!");

}*/