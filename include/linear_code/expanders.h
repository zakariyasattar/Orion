#pragma once

#include "parameter.h"
#include <vector>
#include <span>
#include <iostream>
#include <iomanip>
#include <linear_gkr/prime_field.h>
class graph
{
public:
    int degree;
    std::vector<std::vector<long long>> neighbor;
    std::vector<std::vector<long long>> r_neighbor;
    std::vector<std::vector<prime_field::field_element>> weight;
    std::vector<std::vector<prime_field::field_element>> r_weight;
    long long L, R;
};

extern graph C[100], D[100];

inline void print_adj_mat(std::span<graph> graphs) {
    for(auto& graph : graphs) {
        if(graph.neighbor.size() == 0) continue;
        
        std::cout << "------------- START GRAPH (L=" << graph.L << ", R=" << graph.R << ", degree=" << graph.degree << ") -------------" << std::endl;
        
        std::vector<std::vector<long long>> adj_matrix(graph.L + graph.R,
                                                        std::vector<long long>(graph.L + graph.R, 0));
        
        // Build adjacency matrix
        for(long long i = 0; i < static_cast<long long>(graph.neighbor.size()); i++) {
            std::vector<long long>& neighbors = graph.neighbor[i];
            for(long long j = 0; j < static_cast<long long>(neighbors.size()); j++) {
                adj_matrix[i][graph.L + neighbors[j]] = 1;
            }
        }
        
        // Print column numbers
        std::cout << "     ";  // Space for row labels
        for(long long col = 0; col < graph.L + graph.R; col++) {
            if(col == graph.L) std::cout << "| ";
            std::cout << std::setw(2) << col << " ";
        }
        std::cout << std::endl;
        
        // Print L/R labels for columns
        std::cout << "     ";
        for(long long col = 0; col < graph.L; col++) {
            std::cout << " L ";
        }
        std::cout << "| ";
        for(long long col = 0; col < graph.R; col++) {
            std::cout << " R ";
        }
        std::cout << std::endl;
        
        // Print separator line
        std::cout << "     ";
        for(long long col = 0; col < graph.L; col++) {
            std::cout << "---";
        }
        std::cout << "+-";
        for(long long col = 0; col < graph.R; col++) {
            std::cout << "---";
        }
        std::cout << std::endl;
        
        // Print rows with data
        for(long long node = 0; const auto& row : adj_matrix) {
            if(node < graph.L) {
                std::cout << "L" << std::setw(2) << std::left << node << std::right << ": ";
            } else {
                std::cout << "R" << std::setw(2) << std::left << (node - graph.L) << std::right << ": ";
            }
            
            for(long long col = 0; col < graph.L + graph.R; col++) {
                if(col == graph.L) std::cout << "| ";
                std::cout << std::setw(2) << row[col] << " ";
            }
            std::cout << std::endl;
            node++;
        }
        
        std::cout << "------------- END GRAPH -------------" << std::endl;
    }
}

inline graph generate_random_expander(long long L, long long R, long long d)
{
    graph ret;
    ret.degree = d;
    ret.neighbor.resize(L);
    ret.weight.resize(L);

    ret.r_neighbor.resize(R);
    ret.r_weight.resize(R);
    for(long long i = 0; i < L; ++i)
    {
        ret.neighbor[i].resize(d);
        ret.weight[i].resize(d);
        for(long long j = 0; j < d; ++j)
        {
            long long target = rand() % R;
            prime_field::field_element weight = prime_field::random();
            ret.neighbor[i][j] = target;
            ret.r_neighbor[target].push_back(i);
            ret.r_weight[target].push_back(weight);
            ret.weight[i][j] = weight;
        }
    }
    ret.L = L;
    ret.R = R;
    return ret;
}

inline long long expander_init(long long n, int dep = 0)
{
    //random graph
    if(n <= distance_threshold)
    {
        return n;
    }
    else
    {
        C[dep] = generate_random_expander(n, (long long)(alpha * n), cn);
        long long L = expander_init((long long)(alpha * n), dep + 1);
        D[dep] = generate_random_expander(L, (long long)(n * (r - 1) - L), dn);
        return n + L + (long long)(n * (r - 1) - L);
    }
}