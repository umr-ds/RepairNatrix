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

/*
  a = pointer to the start of the first substring,
  b = pointer to the start of the second substring,
  n = length of the substring to compare
*/
static bool substring_equal(char *a, char *b, int n) { // or use strncmp() == 0 ?
   int i = 0;
   while (i < n) {
       if (a[i] != b[i]) {
           return false;
       }
       i += 1;
   }
   return true;
}


static PyObject* kmer_counting_error_val(PyObject* self, PyObject *args)
{
//TODO: add initial kmer to the result if it IS a violating kmer + inspect why numbers increase so much!
    char* seq;
    unsigned int k, upper_bound = 0;
    if (!PyArg_ParseTuple(args, "sii", &seq, &k, &upper_bound)) {
        return NULL;
    }
    unsigned int seq_len = strlen(seq);
    double *res = calloc(seq_len, sizeof(double));
    double *current_kmer_tmp_result = calloc(seq_len, sizeof(double));
    bool *seen_kmer_start_points = calloc(seq_len, sizeof(bool));
    for (unsigned int i = 0; i <= (seq_len - k); i++) {
        if (seen_kmer_start_points[i]) {
            // this kmer was already seen, we can safely skip it...
            continue;
        }
        unsigned int count = 0;
        for (unsigned int j = i; j <= seq_len - k; j += 1) {
            // compare seq[i:i+k] with seq[j:j+k]
            if (substring_equal(&seq[i], &seq[j], k)) {
                // if they are equal, we set seen to true to skip this position in the future
                seen_kmer_start_points[j] = true;
                // we increase the number of occurrences of the current kmer
                count++;
                for (unsigned int h = 0; h < k; h++) {
                    // increase the kmer-occurrence for each base in the current kmer ( seq[j:j+h] )
                    current_kmer_tmp_result[j + h] += 1;
                }
            }
        }
        // if we got more than upper_bound occurrences of this kmer, we add it to the result
        if (count >= upper_bound) {
            // we got a kmer > upper_bound: update result
            for (unsigned int j = i; j < seq_len; j++) {
                if (current_kmer_tmp_result[j] > 0)
                    res[j] += current_kmer_tmp_result[j] * (count - upper_bound);
            }
        }
        // reset the temporary array to all 0;
        for (unsigned int j = i; j < strlen(seq); j++) {
            current_kmer_tmp_result[j] = 0;
        }
    }
    npy_intp dims[1];
    dims[0] = seq_len;
    PyObject *out = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void *)res);
    free(current_kmer_tmp_result);
    free(seen_kmer_start_points);
    PyArray_ENABLEFLAGS((PyArrayObject*)out, NPY_ARRAY_OWNDATA);
    return out;
}


static PyObject* kmer_counting_error_val_OLD(PyObject* self, PyObject *args)
{
//TODO: add initial kmer to the result if it IS a violating kmer + inspect why numbers increase so much!
  char* seq;
  unsigned int k, upper_bound = 0;
  if (!PyArg_ParseTuple(args, "sii", &seq, &k, &upper_bound)) {
       return NULL;
 }
  unsigned int seq_len = strlen(seq);
  // we store the position of all duplicated substrings in this array. (temporarily even for k mers below the upper bound))
  double* res = calloc(seq_len, sizeof(double));
  // create an array:
  unsigned int* kmer_pos = calloc((seq_len+1), sizeof(int));
  // only if we have more kmers than we have nucleotides, we will have a problem... (wont happen since k >= 1)
  unsigned int kmer_pos_current_size = 0; // keep track of the next free position in the array
  for (unsigned int i = 0; i < (seq_len-k+1); i++) {
    //if i is equal to any of the positions stored in kmer_post, we can skip it
    for (unsigned int j = 0; j < kmer_pos_current_size; j++) {
      if (kmer_pos[j] == i) {
        PySys_WriteStdout("\nSkipping kmer at position %d\n", i);
        goto cnt;
      } else {
        PySys_WriteStdout("%d, ", kmer_pos[j]);
      }
    }
    // backup the current size of the array in case the kmer is below the upper_bound
    unsigned int kmer_pos_current_size_backup = kmer_pos_current_size;
    // slide over the remaining sequence to find other occurences of the current kmer:
    unsigned int count = 0;
    for (unsigned int j = i; j < seq_len-k+1; j++) {
      if (substring_equal(&seq[i], &seq[j], k)) {
        count++;
      }
    }
    for (unsigned int j = i; j < seq_len-k+1; j++) {
      if (substring_equal(&seq[i], &seq[j], k)) {
        PySys_WriteStdout("substr's are equal: j=%d, i=%d, k=%d\n", j, i, k);
        kmer_pos[kmer_pos_current_size] = count;
        kmer_pos_current_size += 1;
      }
    }
    // check how many times the kmer occurs in the sequence:
    unsigned int kmer_size = kmer_pos_current_size - kmer_pos_current_size_backup;
    if (kmer_size > upper_bound) {
      // we found a kmer that is above the upper bound
      for (unsigned int i = kmer_pos_current_size_backup; i < kmer_pos_current_size; i++) {
        for (unsigned int j = 0; j < k; j++) {
          res[i+j] = 1.0 * kmer_size - upper_bound;
        }
      }
    } else {
        //reset the current position in the array...
        kmer_pos_current_size = kmer_pos_current_size_backup;
    }
  cnt:;
  }
  //convert to PyArray:
  npy_intp dims[1];
  dims[0] = seq_len;
  //TODO find out why the first 2 values are 0 instead of the correct values!
  PyObject *out = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, (void *)res);
  free(res);
  free(kmer_pos);
  return out;
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
static char kmer_counting_error_value_docs[] =
"kmer_counting_error_val(seq,k, upper_bound): returns the error value for all kmers of length k in seq with a upper bound of upper_bound";
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
   {"kmer_counting_error_val", kmer_counting_error_val, METH_VARARGS, kmer_counting_error_value_docs},
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