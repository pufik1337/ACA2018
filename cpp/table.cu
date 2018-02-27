/* Copyright (c) 2010-2011, Panos Louridas, GRNET S.A.
 
   All rights reserved.
  
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
 
   * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
 
   * Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the
   distribution.
 
   * Neither the name of GRNET S.A, nor the names of its contributors
   may be used to endorse or promote products derived from this
   software without specific prior written permission.
  
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
   HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
   OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <map>
#include <math.h>
#include <string>
#include <cstring>
#include <limits>
#include "stdio.h"
#include "table.cuh"


void Table::reset() {
    num_outgoing.clear();
    rows.clear();
    nodes_to_idx.clear();
    idx_to_nodes.clear();
    pr.clear();
}

Table::Table(double a, double c, size_t i, bool t, bool n, string d)
    : trace(t),
      alpha(a),
      convergence(c),
      max_iterations(i),
      delim(d),
      numeric(n) {
}

void Table::reserve(size_t size) {
    num_outgoing.reserve(size);
    rows.reserve(size);
}

const size_t Table::get_num_rows() {
    return rows.size();
}

void Table::set_num_rows(size_t num_rows) {
    num_outgoing.resize(num_rows);
    rows.resize(num_rows);
}

const void Table::error(const char *p,const char *p2) {
    cerr << p <<  ' ' << p2 <<  '\n';
    exit(1);
}

const double Table::get_alpha() {
    return alpha;
}

void Table::set_alpha(double a) {
    alpha = a;
}

const unsigned long Table::get_max_iterations() {
    return max_iterations;
}

void Table::set_max_iterations(unsigned long i) {
    max_iterations = i;
}

const double Table::get_convergence() {
    return convergence;
}

void Table::set_convergence(double c) {
    convergence = c;
}

const vector<double>& Table::get_pagerank() {
    return pr;
}

const string Table::get_node_name(size_t index) {
    if (numeric) {
        stringstream s;
        s << index;
        return s.str();
    } else {
        return idx_to_nodes[index];
    }
}

const map<size_t, string>& Table::get_mapping() {
    return idx_to_nodes;
}

const bool Table::get_trace() {
    return trace;
}

void Table::set_trace(bool t) {
    trace = t;
}

const bool Table::get_numeric() {
    return numeric;
}

void Table::set_numeric(bool n) {
    numeric = n;
}

const string Table::get_delim() {
    return delim;
}

void Table::set_delim(string d) {
    delim = d;
}

/*
 * From a blog post at: http://bit.ly/1QQ3hv
 */
void Table::trim(string &str) {

    size_t startpos = str.find_first_not_of(" \t");

    if (string::npos == startpos) {
        str = "";
    } else {
        str = str.substr(startpos, str.find_last_not_of(" \t") - startpos + 1);
    }
}

size_t Table::insert_mapping(const string &key) {

    size_t index = 0;
    map<string, size_t>::const_iterator i = nodes_to_idx.find(key);
    if (i != nodes_to_idx.end()) {
        index = i->second;
    } else {
        index = nodes_to_idx.size();
        nodes_to_idx.insert(pair<string, size_t>(key, index));
        idx_to_nodes.insert(pair<size_t, string>(index, key));;
    }

    return index;
}

int Table::read_file(const string &filename) {

    pair<map<string, size_t>::iterator, bool> ret;

    reset();
    
    istream *infile;

    if (filename.empty()) {
      infile = &cin;
    } else {
      infile = new ifstream(filename.c_str());
      if (!infile) {
          error("Cannot open file", filename.c_str());
      }
    }
    
    size_t delim_len = delim.length();
    size_t linenum = 0;
    string line; // current line
    while (getline(*infile, line)) {
        string from, to; // from and to fields
        size_t from_idx, to_idx; // indices of from and to nodes
        size_t pos = line.find(delim);
        if (pos != string::npos) {
            from = line.substr(0, pos);
            trim(from);
            if (!numeric) {
                from_idx = insert_mapping(from);
            } else {
                from_idx = strtol(from.c_str(), NULL, 10);
            }
            to = line.substr(pos + delim_len);
            trim(to);
            if (!numeric) {
                to_idx = insert_mapping(to);
            } else {
                to_idx = strtol(to.c_str(), NULL, 10);
            }
            add_arc(from_idx, to_idx);
        }

        linenum++;
        if (linenum && ((linenum % 100000) == 0)) {
            cerr << "read " << linenum << " lines, "
                 << rows.size() << " vertices" << endl;
        }

        from.clear();
        to.clear();
        line.clear();
    }

    cerr << "read " << linenum << " lines, "
         << rows.size() << " vertices" << endl;

    nodes_to_idx.clear();

    if (infile != &cin) {
        delete infile;
    }
    reserve(idx_to_nodes.size());

    //vector<unsigned> row_rep = {0};
    //vector<unsigned> col_rep;
    row_rep.push_back(0);

    for (unsigned nodes = 0; nodes < rows.size(); nodes++){
        acc += rows[nodes].size();
        row_rep.push_back(acc);
        col_rep.insert(col_rep.end(), rows[nodes].begin(), rows[nodes].end());
        //if (rows[nodes].size() == 0)
            //zero_outgoing_idxvec.push_back(nodes);
    }

    return 0;
}

/*
 * Taken from: M. H. Austern, "Why You Shouldn't Use set - and What You Should
 * Use Instead", C++ Report 12:4, April 2000.
 */
template <class Vector, class T>
bool Table::insert_into_vector(Vector& v, const T& t) {
    typename Vector::iterator i = lower_bound(v.begin(), v.end(), t);
    if (i == v.end() || t < *i) {
        v.insert(i, t);
        return true;
    } else {
        return false;
    }
}

bool Table::add_arc(size_t from, size_t to) {

    bool ret = false;
    size_t max_dim = max(from, to);
    if (trace) {
        cout << "checking to add " << from << " => " << to << endl;
    }
    if (rows.size() <= max_dim) {
        max_dim = max_dim + 1;
        if (trace) {
            cout << "resizing rows from " << rows.size() << " to "
                 << max_dim << endl;
        }
        rows.resize(max_dim);
        if (num_outgoing.size() <= max_dim) {
            num_outgoing.resize(max_dim);
        }
    }

    ret = insert_into_vector(rows[to], from);

    if (ret) {
        num_outgoing[from]++;
        if (trace) {
            cout << "added " << from << " => " << to << endl;
        }
    }

    return ret;
}

__global__ void danglingPr();

__global__ void pageRank(
    const size_t* csrRowPtrA, 
    const size_t* csrColPtrA, 
    const size_t* num_outgoing, 
    double* pr, 
    double* old_pr,
    double* diff, 
    const double alpha,
    const double one_Av,
    const double one_Iv,
    const double sum_pr)
    {
    
        unsigned ki = blockIdx.x*blockDim.x + threadIdx.x;
        old_pr[ki] = pr[ki]/sum_pr;
        unsigned row_start = csrRowPtrA[ki];
        unsigned row_end = csrRowPtrA[ki+1];
        //printf ("Thread number %d; Row Start: %d; Row end = %d \n", threadIdx.x, row_start, row_end);
        double h = 0.0;
        for(unsigned i = row_start; i < row_end; i++){
            unsigned ci = csrColPtrA[i];
            double h_v = (num_outgoing[ci])
            ? 1.0 / num_outgoing[ci]
            : 0.0;
            h += h_v * old_pr[ci];
            //printf ("Thread number %d; Column index (pointing at me): %d; h = %f \n", threadIdx.x, ci, h);       
        }
        h *= alpha;
        double this_pr = h + one_Av + one_Iv; 
        pr[ki] = this_pr;
        diff[ki] = fabs(this_pr - old_pr[ki]);
        //pr[ki] = ki;
}

void Table::pagerank() {

    vector<size_t>::iterator ci; // current incoming
    size_t num_rows = rows.size();

    double diff = 1;
    size_t i;
    double sum_pr; // sum of current pagerank vector elements
    double dangling_pr = 0; // sum of current pagerank vector elements for dangling
    			// nodes
    unsigned long num_iterations = 0;
    vector<double> old_pr;



    //print_table();

    //print_outgoing();

    
    if (num_rows == 0) {
        return;
    }
    
    pr.resize(num_rows);

    pr[0] = 1;

    if (trace) {
        print_pagerank();
    }
  
    vector<unsigned> zero_outgoing_idxvec;

    for (size_t k = 0; k < pr.size(); k++) {
        if (num_outgoing[k] == 0) {
            zero_outgoing_idxvec.push_back(k);
        }
    }

    thrust::device_vector<unsigned> d_zero_outgoing_index;
    d_zero_outgoing_index = zero_outgoing_idxvec;

    // allocate num_outgoing, pr, old_pr, csrRowPtrA, csrColPtrA
    size_t int_size = num_rows * sizeof(size_t);
    size_t double_size = num_rows * sizeof(double);
    size_t col_size = col_rep.size() * sizeof(size_t);
    size_t row_size = (row_rep.size() + 1) * sizeof(size_t);
    size_t dangling_size = (zero_outgoing_idxvec.size()) * sizeof(double);

    size_t* d_num_outgoing;
    size_t* d_row_rep;
    size_t* d_col_rep;
    double* d_pr;
    double* d_old_pr;
    double* d_diff;
    double* d_dangling;

    cudaMalloc(&d_num_outgoing, int_size);
    cudaMalloc(&d_row_rep, row_size);
    cudaMalloc(&d_col_rep, col_size);
    cudaMalloc(&d_pr, double_size);
    cudaMalloc(&d_old_pr, double_size);
    cudaMalloc(&d_diff, double_size);
    cudaMalloc(&d_dangling, dangling_size);

    cudaMemcpy(d_num_outgoing, &num_outgoing[0], int_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_row_rep, &row_rep[0], row_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_col_rep, &col_rep[0], col_size, cudaMemcpyHostToDevice);
    //cudaMemcpy(d_old_pr, &pr[0], double_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pr, &pr[0], double_size, cudaMemcpyHostToDevice);

    double one_Av;
    double one_Iv = (1 - alpha) / num_rows;

    while (diff > convergence && num_iterations < max_iterations) {
        
        thrust::gather(thrust::device,
            d_zero_outgoing_index.begin(), d_zero_outgoing_index.end(),
            d_pr,
            d_dangling
        );

        dangling_pr = thrust::reduce(thrust::device, d_dangling, d_dangling + zero_outgoing_idxvec.size());
        sum_pr = thrust::reduce(thrust::device, d_pr, d_pr + num_rows); 
                             
        one_Av = alpha * dangling_pr / num_rows;

        pageRank<<<38000, 500>>>(
            d_row_rep,
            d_col_rep, 
            d_num_outgoing,
            d_pr,
            d_old_pr,
            d_diff,
            alpha,
            one_Av,
            one_Iv,
            sum_pr);

        cudaDeviceSynchronize();

        diff = thrust::reduce(thrust::device, d_diff, d_diff + num_rows);

        num_iterations++;

    }

    cudaMemcpy(&pr[0], d_pr, double_size, cudaMemcpyDeviceToHost);

    cudaFree(d_col_rep);
    cudaFree(d_row_rep);
    cudaFree(d_num_outgoing);
    cudaFree(d_old_pr);
    cudaFree(d_pr);
    cudaFree(d_dangling);
}

const void Table::print_params(ostream& out) {
    out << "alpha = " << alpha << " convergence = " << convergence
        << " max_iterations = " << max_iterations
        << " numeric = " << numeric
        << " delimiter = '" << delim << "'" << endl;
}

const void Table::print_table() {
    vector< vector<size_t> >::iterator cr;
    vector<size_t>::iterator cc; // current column

    size_t i = 0;
    for (cr = rows.begin(); cr != rows.end(); cr++) {
        cout << i << ":[ ";
        for (cc = cr->begin(); cc != cr->end(); cc++) {
            if (numeric) {
                cout << *cc << " ";
            } else {
                cout << idx_to_nodes[*cc] << " ";
            }
        }
        cout << "]" << endl;
        i++;
    }
}

const void Table::print_outgoing() {
    vector<size_t>::iterator cn;

    cout << "[ ";
    for (cn = num_outgoing.begin(); cn != num_outgoing.end(); cn++) {
        cout << *cn << " ";
    }
    cout << "]" << endl;

}

const void Table::print_row_col() {
    vector<size_t>::iterator cn;

    cout << "Row: \n [ ";
    for (cn = row_rep.begin(); cn != row_rep.end(); cn++) {
        cout << *cn << " ";
    }
    cout << "]" << endl;

    cout << "Column: \n [ ";
    for (cn = col_rep.begin(); cn != col_rep.end(); cn++) {
        cout << *cn << " ";
    }
    cout << "]" << endl;

}

const void Table::print_pagerank() {

    vector<double>::iterator cr;
    double sum = 0;

    cout.precision(numeric_limits<double>::digits10);
    
    cout << "(" << pr.size() << ") " << "[ ";
    for (cr = pr.begin(); cr != pr.end(); cr++) {
        cout << *cr << " ";
        sum += *cr;
        cout << "s = " << sum << " ";
    }
    cout << "] "<< sum << endl;
}

const void Table::print_pagerank_v() {

    size_t i;
    size_t num_rows = pr.size();
    double sum = 0;
    
    cout.precision(numeric_limits<double>::digits10);

    for (i = 0; i < num_rows; i++) {
        if (!numeric) {
            cout << idx_to_nodes[i] << " = " << pr[i] << endl;
        } else {
            cout << i << " = " << pr[i] << endl;
        }
        sum += pr[i];
    }
    cerr << "s = " << sum << " " << endl;
}
