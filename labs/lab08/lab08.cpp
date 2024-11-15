#include <iostream>
#include <fstream>
#include <string>
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
#include <algorithm>    // std::min_element, std::max_element
#include <array>
#include <mpi.h>

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

  string info(const string name="") const
  {
    return (name.empty() ? "CSR" : name) + "<" + to_string(nrows) + ", " + to_string(ncols) + ", " +
      (ptr ? to_string(ptr[nrows]) : "0") + ">";
  }

  /** 
   * Read the matrix from a CLUTO file.
   * The first line is "nrows ncols nnz".
   * Each other line is a row in the matrix containing column ID-value pairs for non-zeros in the row.
  */
  void read(const std::string &filename)
  {
    FILE * infile = fopen(filename.c_str(), "r");
    char * line = NULL;
    size_t n, nr, nnz;
    char *head;
    char *tail;
    idx_t cid;
    double dval;
    
    if (!infile) {
      throw std::runtime_error("Could not open CLU file\n");
    }
    if(getline (&line, &n, infile) < 0){
      throw std::runtime_error("Could not read first line from CLU file\n");
    }
    //read matriz size info
    size_t rnrows, rncols, rnnz;
    sscanf(line, "%zu %zu %zu", &rnrows, &rncols, &rnnz);

    //allocate space
    this->reserve(rnrows, rnnz);
    ncols = rncols;
    
    //read in rowval, rowind, rowptr
    this->ptr[0]= 0;
    nnz = 0;
    nr = 0;

    while(getline(&line, &n, infile) != -1){
      head = line;
      while (1) {
        cid = (idx_t) strtol(head, &tail, 0);
        if (tail == head)
          break;
        head = tail;

        if(cid <= 0){
          throw std::runtime_error("Invalid column ID while reading CLUTO matrix\n");
        }
        this->ind[nnz] = cid - 1; //csr/clu files have 1-index based column IDs and our matrix is 0-based.
        dval = strtod(head, &tail);
        head = tail;
        this->val[nnz++] = dval;
      }
      this->ptr[nr+1] = nnz;
      nr++;
    }
    assert(nr == rnrows);
    free(line);
    fclose(infile);
  }
  
  /** 
   * Read the matrix from a CLUTO file.
   * The first line is "nrows ncols nnz".
   * Each other line is a row in the matrix containing column ID-value pairs for non-zeros in the row.
  */
  static csr_t * from_CLUTO(const std::string &filename)
  {
    auto mat = new csr_t();
    mat->read(filename);
    return mat;
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

  /**
   * Normalize the rows of the matrix by L1 or L2 norm.
   * @param norm The norm of the matrix to normalize by.
  */
  void normalize(int norm=2)
  {
    #pragma omp parallel if (ptr[nrows] > 1e+6)
    {
      val_t sum;
      #pragma omp for schedule (static)
      for (idx_t i = 0; i < nrows; i++) { // each row
        sum = 0;
        for (ptr_t j = ptr[i]; j < ptr[i + 1]; j++) { // each value in row
          if (norm == 2) {
            sum += val[j] * val[j];
          } else if (norm == 1) {
            sum += val[j] > 0 ? val[j] : -val[j];
          } else {
            throw std::runtime_error("Norm must be 1 or 2.");
          }
        }
        if (sum > 0) {
          if (norm == 2) {
            sum = (double) 1.0 / sqrt(sum);
          } else {
            sum = (double) 1.0 / sum;
          }
          for (ptr_t j = ptr[i]; j < ptr[i + 1]; j++) {
            val[j] *= sum;
          }
        }
      }
    }
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
 * Read a subset of the rows in the CSR matrix. 
 * Stop reading at the end of the row once at least rnnz non-zeros have been read.
 * Note that the file pointer may be in the middle of the matrix or right after reading the "nrows ncols nnz" line.
 * 
 * @param infile Pointer to the opened input file
 * @param rnnz   Number of non-zeros to read
 * @param nrows  Expected number of rows (may be more or less)
 * @param ncols  Number of columns in the matrix
*/
csr_t * read_partial_csr(FILE * infile, ptr_t rnnz, idx_t nrows, idx_t ncols)
{
  // TO IMPLEMENT




  // don't forget to reduce the space of the matrix to only what is needed

}

/** 
 * Read the first line in a cluto file and return nrows, ncols, and nnz.
 * 
 * @param infile path to the input file
 * @param nrows  Pointer to where nrows should be stored
 * @param ncols  Pointer to where ncols should be stored
 * @param nnz    Pointer to where nnz should be stored
*/
void get_cluto_stats(FILE * infile, idx_t * nrows, idx_t * ncols, ptr_t * nnz)
{
  char * line = NULL;
  size_t n;
  
  if (!infile) {
    throw std::runtime_error("Invalid file pointer\n");
  }
  if(getline (&line, &n, infile) < 0){
    throw std::runtime_error("Could not read first line from CLU file\n");
  }
  //read matriz size info
  sscanf(line, "%u %u %zu", nrows, ncols, nnz);
  free(line);
}

/** 
 * Read the first line in a cluto file and return nrows, ncols, and nnz.
 * 
 * @param infile Pointer to the opened input file
 * @param nrows  Pointer to where nrows should be stored
 * @param ncols  Pointer to where ncols should be stored
 * @param nnz    Pointer to where nnz should be stored
*/
void get_cluto_stats(char * fname, idx_t * nrows, idx_t * ncols, ptr_t * nnz)
{
  FILE * infile = fopen(fname, "r");
  get_cluto_stats(infile, nrows, ncols, nnz);
  fclose(infile);
}



/**
 * Send the csr matrix mat to process rank.
 * @param mat CSR matrix to send
 * @param rank Process ID to send to
 * @param comm Communicator to send via
 */
void send_csr(csr_t * mat, int to, MPI_Comm comm) {
  // TO IMPLEMENT
}

/**
 * Receive a csr matrix from process rank.
 * @param rank Process ID receiving matrix from
 * @param comm Communicator to receive via
 */
csr_t * receive_csr(int rank, MPI_Comm comm) {
  
  // TO IMPLEMENT

  // return the received matrix
  return mat;
}


/**
 * Read a CLUTO matrix from a file and split it across all the processes in the communicator.
 */
csr_t * read_and_bcast_csr(char * fname, MPI_Comm comm)
{
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(comm, &world_rank);
  int world_size;
  MPI_Comm_size(comm, &world_size);

  // TO IMPLEMET

  // return the local subset of the matrix for the current process
  return mat;
}

int main(int argc, char *argv[])
{
 
  if(argc < 2){
    cerr << "Invalid options." << endl << "<program> <cluto_matrix>" << endl;
    exit(1);
  }

  MPI_Init(&argc, &argv);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  // read and distribute the input cluto matrix
  if(world_rank == 0) cout << "Process 0: Broadcasting DB matrix." << endl;
  auto mat = read_and_bcast_csr(argv[1], MPI_COMM_WORLD); // local subset of DB matrix
  cout << "Process " << world_rank << " Matrix subset info: " << mat->info() << endl;

  // make sure all rows and non-zeros were distributed correclty
  idx_t nrows, ncols, tnrows;
  ptr_t nnz, tnnz;
  get_cluto_stats(argv[1], &nrows, &ncols, &nnz);
  MPI_Allreduce(&mat->nrows, &tnrows, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&mat->ptr[mat->nrows], &tnnz, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);
  if(world_rank == 0)
    if(tnrows != nrows || tnnz != nnz){
      cerr << "Error: Distributed matrix has " << tnrows << " rows and " << tnnz << " non-zeros, but should have " << nrows << " rows and " << nnz << " non-zeros." << endl;
    } else {
      cout << "Distribution successful!" << endl;
      cout << "Global matrix stats: CSR<" << nrows << ", " << ncols << ", " << nnz << ">" << endl; 
    }

  delete mat;

  // Finalize the MPI environment.
  MPI_Finalize();

  return 0;
}

