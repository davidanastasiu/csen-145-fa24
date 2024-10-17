#include <fstream>
#include <iostream>
#include <omp.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cstring>      /* strcasecmp */
#include <cstdint>
#include <assert.h>
#include <vector>       // std::vector
#include <algorithm>    // std::random_shuffle
#include <random>
#include <stdexcept>

using namespace std;

using idx_t = std::uint32_t;
using val_t = float;
using ptr_t = std::uintptr_t;

/**
 * CSR structure to store search results
 */
typedef struct csr_t {
  idx_t nrows; // number of rows
  idx_t ncols; // number of rows
  idx_t * ind; // column ids
  val_t * val; // values
  ptr_t * ptr; // pointers (start of row in ind/val)

  csr_t()
  {
    nrows = ncols = 0;
    ind = nullptr;
    val = nullptr;
    ptr = nullptr;
  }

  /**
   * Reserve space for more rows or non-zeros. Structure may only grow, not shrink.
   * @param nrows Number of rows
   * @param nnz   Number of non-zeros
   */
  void reserve(const idx_t nrows, const ptr_t nnz)
  {
    if(nrows > this->nrows){
      if(ptr){
        ptr = (ptr_t*) realloc(ptr, sizeof(ptr_t) * (nrows+1));
      } else {
        ptr = (ptr_t*) malloc(sizeof(ptr_t) * (nrows+1));
        ptr[0] = 0;
      }
      if(!ptr){
        throw std::runtime_error("Could not allocate ptr array.");
      }
    }
    if(ind){
      ind = (idx_t*) realloc(ind, sizeof(idx_t) * nnz);
    } else {
      ind = (idx_t*) malloc(sizeof(idx_t) * nnz);
    }
    if(!ind){
      throw std::runtime_error("Could not allocate ind array.");
    }
    if(val){
      val = (val_t*) realloc(val, sizeof(val_t) * nnz);
    } else {
      val = (val_t*) malloc(sizeof(val_t) * nnz);
    }
    if(!val){
      throw std::runtime_error("Could not allocate val array.");
    }
    this->nrows = nrows;
  }

  csr_t ( const csr_t &other)
  {
    this->nrows = this->ncols = 0;
    this->ptr = nullptr;
    this->ind = nullptr;
    this->val = nullptr;
    this->reserve(other.nrows, other.ptr[other.nrows]);
    memcpy(ptr, other.ptr, sizeof(ptr_t) * (nrows+1));
    memcpy(ind, other.ind, sizeof(idx_t) * ptr[nrows]);
    memcpy(val, other.val, sizeof(val_t) * ptr[nrows]);
    this->ncols = other.ncols;
  }

  /** 
   * Transpose matrix
   */
  void transpose()
  {
    /* TODO: Implement matrix transposition */
  }

  /**
   * Create random matrix with given sparsity factor.
   * @param nrows Number of rows
   * @param ncols Number of columns
   * @param factor   Sparsity factor
   */
  static csr_t * random(const idx_t nrows, const idx_t ncols, const double factor)
  {
    ptr_t nnz = (ptr_t) (factor * nrows * ncols);
    if(nnz >= nrows * ncols / 2.0){
      throw std::runtime_error("Asking for too many non-zeros. Matrix is not sparse.");
    }
    auto mat = new csr_t();
    mat->reserve(nrows, nnz);
    mat->ncols = ncols;

    // create a vector of random numbers of size nrows
    vector<double> numbers(nrows);
    std::random_device r;
    std::default_random_engine gen(r());
    std::uniform_int_distribution<int> uniform_dist(1, ncols * factor);
    std::generate(numbers.begin(), numbers.end(), [&] {return uniform_dist(gen); });
    auto sum = std::accumulate(numbers.begin(), numbers.end(), 0.0);
    // normalize the vector to sum to 1
    for (idx_t i = 0; i < nrows; i++) {
      numbers[i] = numbers[i] / (double) sum;
    }
    /* fill in ptr array; generate random row sizes */
    for(idx_t i=0; i < mat->nrows; ++i){
      mat->ptr[i+1] = mat->ptr[i] + std::max((ptr_t)(numbers[i] * nnz), (ptr_t)1);
      if(mat->ptr[i+1] > nnz){
        mat->ptr[i+1] = nnz;
      }
    }

    /* fill in indices and values with random numbers */
    #pragma omp parallel
    {
      int tid = omp_get_thread_num();
      unsigned int seed = (unsigned long) mat * (1+tid);
      std::vector<int> perm;
      for(idx_t i=0; i < ncols; ++i){
        perm.push_back(i);
      }
      std::random_device seeder;
      std::mt19937 engine(seeder());

      #pragma omp for
      for(idx_t i=0; i < nrows; ++i){
        std::shuffle(perm.begin(), perm.end(), engine);
        for(ptr_t j=mat->ptr[i]; j < mat->ptr[i+1]; ++j){
          mat->ind[j] = perm[j - mat->ptr[i]];
          mat->val[j] = ((double) rand_r(&seed)/rand_r(&seed));
        }
      }
    }

    return mat;
  }

  /* check if matrix is equal to the other matrix */
  bool equal(csr_t * other, const double eps = 1e-3)
  {
    if(nrows != other->nrows || ncols != other->ncols || ptr[nrows] != other->ptr[other->nrows]){
      return false;
    }
    // matrices need to be sorted to check equality
    if(!this->indices_sorted()){
      this->sort_indices();
    }
    if(!other->indices_sorted()){
      other->sort_indices();
    }
    for(idx_t i=0; i < nrows; ++i){
      for(ptr_t j=ptr[i]; j < ptr[i+1]; ++j){
        if(ind[j] != other->ind[j]){
          cout << "row " << i << " cid missmatch: " << ind[j] << " != " << other->ind[j] << endl;
          return false;
        }
        if(abs(val[j] - other->val[j]) > eps){
          cout << "row " << i << " val missmatch: " << val[j] << " != " << other->val[j] << endl;
          return false;
        }
      }
    }
    return true;
  }

  /* sort non-zeros in each row in increasing order of columnd ids */
  void sort_indices()
  {
    /* TODO: Implement index sorting */
  }

  /* check if indices are sorted */
  bool indices_sorted()
  {
    for(idx_t i=0; i < nrows; ++i){
      for(ptr_t j=ptr[i]; j < ptr[i+1]-1; ++j){
        if(ind[j+1] <= ind[j]){
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Write matrix to text file
   * @param output_fpath File to write to
   */
  void write(const std::string output_fpath, const bool header=false)
  {
    std::fstream resfile;
    resfile.open(output_fpath, std::ios::out);
    if(!resfile){
      throw std::runtime_error("Could not open output file for writing.");
    }
    if(header){
      resfile << nrows << " " << ncols << " " << ptr[nrows] << std::endl;
    }
    for(idx_t i=0; i < nrows; ++i){
      for(ptr_t j=ptr[i]; j < ptr[i+1]; ++j){
        resfile << ind[j] << " " << val[j];
        if(j+1 < ptr[i+1]){
          resfile << " ";
        }
      }
      resfile << std::endl;
    }
    resfile.close();
  }

  string info(const string name="") const
  {
    return (name.empty() ? "CSR" : name) + "<" + to_string(nrows) + ", " + to_string(ncols) + ", " +
      (ptr ? to_string(ptr[nrows]) : "0") + ">";
  }

  ~csr_t()
  {
    if(ind){
      free(ind);
    }
    if(val){
      free(val);
    }
    if(ptr){
      free(ptr);
    }
  }
} csr_t;

/**
 * Ensure the matrix is valid
 * @param mat Matrix to test
 */
void test_matrix(csr_t * mat){
  auto nrows = mat->nrows;
  auto ncols = mat->ncols;
  assert(mat->ptr);
  auto nnz = mat->ptr[nrows];
  for(idx_t i=0; i < nrows; ++i){
    assert(mat->ptr[i] <= nnz);
  }
  for(ptr_t j=0; j < nnz; ++j){
    assert(mat->ind[j] < ncols);
  }
}


/**
 * Multiply A and B and write output in C.
 * Note that C has no data allocations (i.e., ptr, ind, and val pointers are null).
 * Use `csr_t::reserve` to increase C's allocations as necessary.
 * @param A  Matrix A.
 * @param B  Matrix B.
 * @param C  Output matrix
 */
void sparsematmult_sparse_sparse_serial(csr_t * A, csr_t * B, csr_t *C)
{
  // To implement
}

/**
 * Multiply A and B and write output in C.
 * Note that C has no data allocations (i.e., ptr, ind, and val pointers are null).
 * Use `csr_t::reserve` to increase C's allocations as necessary.
 * @param A  Matrix A.
 * @param B  Matrix B.
 * @param C  Output matrix
 */
void sparsematmult_sparse_sparse_parallel(csr_t * A, csr_t * B, csr_t *C)
{
  // To implement
}



int main(int argc, char *argv[])
{
  if(argc < 4){
    cerr << "Invalid options." << endl << "<program> <A_nrows> <A_ncols> <B_ncols> <fill_factor> [-t <num_threads>] [-s]" << endl;
    exit(1);
  }
  int nrows = atoi(argv[1]);
  int ncols = atoi(argv[2]);
  int ncols2 = atoi(argv[3]);
  double factor = atof(argv[4]);
  bool serial = false;
  int nthreads = 1;
  for(int i=5; i < argc; ++i){
    if(strcasecmp(argv[i], "-s") == 0){
      serial = true;
    } else if(strcasecmp(argv[i], "-t") == 0 && i+1 < argc){
      nthreads = atoi(argv[i+1]);
    }
  }
  cout << "A_nrows: " << nrows << endl;
  cout << "A_ncols: " << ncols << endl;
  cout << "B_nrows: " << ncols << endl;
  cout << "B_ncols: " << ncols2 << endl;
  cout << "factor: " << factor << endl;
  cout << "nthreads: " << nthreads << endl;

  /* initialize random seed: */
  srand (time(NULL));

  auto A = csr_t::random(nrows, ncols, factor);
  auto B = csr_t::random(ncols, ncols2, factor); // Note B is not transposed.
  test_matrix(A);
  test_matrix(B);

  cout << A->info("A") << endl;
  cout << B->info("B") << endl;

  auto C = new csr_t(); // Note that C has no data allocations so far.
  if(serial){
    omp_set_num_threads(1);
    auto t1 = omp_get_wtime();
    sparsematmult_sparse_sparse_serial(A, B, C);
    auto t2 = omp_get_wtime();
    cout << C->info("C") << endl;
    cout << "Execution time (serial sparse-sparse dot-product): " << (t2-t1) << endl;
  } else {
    omp_set_num_threads(nthreads);
    omp_set_num_threads(1);
    auto t1 = omp_get_wtime();
    sparsematmult_sparse_sparse_parallel(A, B, C);
    auto t2 = omp_get_wtime();
    cout << C->info("C") << endl;
    cout << "Execution time (parallel sparse-sparse dot-product): " << (t2-t1) << endl;

  }

  delete A;
  delete B;
  delete C;
  
  return 0;
}
