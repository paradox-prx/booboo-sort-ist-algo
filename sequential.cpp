#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <chrono>
#include <stdexcept>

using namespace std;

// Function to compute factorial
long long factorial(int n) {
    long long fact = 1;
    for (int i = 1; i <= n; ++i) fact *= i;
    return fact;
}

// Convert permutation vector to string
string permToString(const vector<int>& perm) {
    string s;
    for (int x : perm) s += to_string(x);
    return s;
}

// Compute r(v): position of the first symbol from the right not in its correct position
int computeR(const vector<int>& v) {
    int n = v.size();
    for (int i = n - 1; i >= 0; i--) {
        if (v[i] != i + 1) return i + 1; // 1-based indexing
    }
    return 0; // Identity permutation
}

// Compute inverse permutation: v^{-1}(i) is the position of symbol i
vector<int> computeInverse(const vector<int>& v) {
    int n = v.size();
    vector<int> inv(n + 1); // 1-based indexing
    for (int i = 0; i < n; i++) {
        inv[v[i]] = i + 1; // Position of symbol v[i]
    }
    return inv;
}

// Swap function: Returns permutation obtained by swapping symbol x with its adjacent symbol
vector<int> Swap(const vector<int>& v, int x) {
    int n = v.size();
    vector<int> p = v;
    int i = computeInverse(v)[x]; // Position of symbol x (1-based)
    if (i >= 1 && i <= n - 1) {
        swap(p[i - 1], p[i]); // Swap symbols at positions i and i+1 (0-based)
    }
    return p;
}

// FindPosition function: Determines the swap position
vector<int> FindPosition(const vector<int>& v, int t) {
    int n = v.size();
    vector<int> p;
    vector<int> identity(n);
    for (int i = 0; i < n; i++) identity[i] = i + 1;

    // Condition 1: t == 2 and Swap(v, t) == 1_n
    p = Swap(v, t);
    if (t == 2 && p == identity) {
        return Swap(v, t - 1);
    }

    // Condition 2: v_{n-1} in {t, n-1}
    if (v[n - 2] == t || v[n - 2] == n - 1) {
        int j = computeR(v);
        return Swap(v, j);
    }

    // Condition 3: Otherwise
    return Swap(v, t);
}

// Parent1 function: Computes parent of vertex v in tree T_t^n
vector<int> Parent1(const vector<int>& v, int t, int n) {
    vector<int> identity(n);
    for (int i = 0; i < n; i++) identity[i] = i + 1;

    // If v is the root (1_n), return empty
    if (v == identity) return {};

    vector<int> p;

    // Case A: v_n == n
    if (v[n - 1] == n) {
        if (t != n - 1) {
            p = FindPosition(v, t);
        } else {
            p = Swap(v, v[n - 2]); // Swap with v_{n-1}
        }
    }
    // Case B: v_n == n-1
    else if (v[n - 1] == n - 1) {
        p = Swap(v, n);
        if (v[n - 2] == n && p != identity) {
            if (t == 1) {
                p = Swap(v, n);
            } else {
                p = Swap(v, t - 1);
            }
        } else {
            if (v[n - 1] == t) {
                p = Swap(v, n);
            } else {
                p = Swap(v, t);
            }
        }
    }
    // Case C: v_n == j for j in {1, 2, ..., n-2}
    else {
        if (v[n - 1] == t) {
            p = Swap(v, n);
        } else {
            p = Swap(v, t);
        }
    }

    return p;
}

// Generate permutations iteratively using std::next_permutation
vector<vector<int>> generatePermutations(int n) {
    vector<vector<int>> perms;
    vector<int> curr(n);
    for (int i = 0; i < n; ++i) curr[i] = i + 1; // Start with 123...n
    do {
        perms.push_back(curr);
    } while (next_permutation(curr.begin(), curr.end()));
    sort(perms.begin(), perms.end()); // Ensure consistent order
    return perms;
}

// Generate the bubble-sort network graph for METIS
void generateGraph(int n, const vector<vector<int>>& perms, vector<idx_t>& xadj, vector<idx_t>& adjncy) {
    xadj.push_back(0);
    for (size_t i = 0; i < perms.size(); ++i) {
        const auto& v = perms[i];
        for (int pos = 1; pos <= n - 1; ++pos) {
            vector<int> neighbor = Swap(v, pos);
            // Find neighbor in perms (binary search since perms is sorted)
            auto it = lower_bound(perms.begin(), perms.end(), neighbor);
            if (it != perms.end() && *it == neighbor) {
                adjncy.push_back(distance(perms.begin(), it));
            }
        }
        xadj.push_back(adjncy.size());
    }
}

// Partition the graph using METIS
void partitionGraph(int n, int nparts, const vector<vector<int>>& perms, vector<idx_t>& part) {
    vector<idx_t> xadj, adjncy;
    generateGraph(n, perms, xadj, adjncy);
    idx_t nvtxs = perms.size();
    idx_t ncon = 1;
    idx_t objval;
    part.resize(nvtxs);
    int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), NULL, NULL, NULL,
                                  &nparts, NULL, NULL, NULL, &objval, part.data());
    if (ret != METIS_OK) {
        throw runtime_error("METIS_PartGraphKway failed");
    }
}

// Get permutations assigned to this rank
vector<vector<int>> getAssignedPermutations(const vector<vector<int>>& perms, int rank, const vector<idx_t>& part) {
    vector<vector<int>> assigned;
    for (size_t i = 0; i < perms.size(); ++i) {
        if (part[i] == rank) {
            assigned.push_back(perms[i]);
        }
    }
    return assigned;
}

// Construct ISTs sequentially (no parallelism)
void constructISTsSequential(int n) {
    if (n > 13) {
        cerr << "Error: n > 13 is not supported due to memory constraints\n";
        return;
    }

    // Generate permutations once
    vector<vector<int>> perms;
    try {
        perms = generatePermutations(n);
    } catch (const bad_alloc& e) {
        cerr << "Error: Memory allocation failed during permutation generation\n";
        return;
    }

    // Partition the graph using METIS
    vector<idx_t> part;
    try {
        partitionGraph(n, 1, perms, part); // Sequential version, no MPI
    } catch (const runtime_error& e) {
        cerr << "Error: " << e.what() << "\n";
        return;
    }

    // Get assigned permutations
    vector<vector<int>> permutations = getAssignedPermutations(perms, 0, part); // Only rank 0 in sequential

    // Open output file for summary
    ofstream out("output_sequential.txt");
    if (!out.is_open()) {
        cerr << "Error: Failed to open output file for sequential run\n";
        return;
    }
    out << "Sequential run: " << permutations.size() << " permutations\n";

    // Measure execution time
    auto start_time = chrono::high_resolution_clock::now();

    // Store results for each tree (minimal output to reduce I/O)
    vector<size_t> pair_counts(n, 0); // Count pairs per tree

    for (int t = 1; t <= n - 1; t++) {
        for (size_t i = 0; i < permutations.size(); i++) {
            vector<int> parent = Parent1(permutations[i], t, n);
            if (!parent.empty()) {
                pair_counts[t]++;
            }
        }
    }

    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end_time - start_time;
    out << "Execution time: " << duration.count() << " seconds\n";

    // Write summary to file
    for (int t = 1; t <= n - 1; t++) {
        out << "Tree T_" << t << "^" << n << ": " << pair_counts[t] << " vertex-parent pairs\n";
    }
    out.close();

    cout << "Results written to output_sequential.txt file\n";
}

int main() {
    int n = 11; // For example, you can change this
    constructISTsSequential(n);

    return 0;
}
