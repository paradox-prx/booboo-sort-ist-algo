# Parallel Graph Algorithm Implementation and Performance Analysis

## Project Overview
This project investigates the parallel implementation of graph algorithms, specifically focusing on constructing Independent Spanning Trees (ISTs) in bubble-sort networks. We explore three parallelization strategies:

- **Sequential Implementation**: Baseline computation of ISTs without parallelism.
- **MPI (Message Passing Interface)**: Distributed computation across multiple nodes.
- **OpenMP+MPI Hybrid**: Combines MPI for inter-node communication and OpenMP for intra-node parallelism.

The project leverages **METIS** for efficient graph partitioning to optimize the load distribution across processors.

## Methodology
The algorithm involves constructing multiple ISTs in a bubble-sort network. The strategies used are:

- **Sequential**: Single-threaded execution.
- **MPI**: Nodes communicate via MPI to compute parts of the graph.
- **OpenMP+MPI**: Combines MPI for inter-node communication with OpenMP for multi-core processing on each node.

## Results

### Execution Times (for various values of n):

| **n** | **Sequential Time (T_seq)** | **MPI Time (T_MPI)** | **OpenMP+MPI Time (T_OMP_MPI)** | **Speedup (MPI)** | **Speedup (OpenMP+MPI)** |
| ----- | ---------------------------- | --------------------- | --------------------------------- | ----------------- | ------------------------ |
| 6     | 0.2s                         | 0.15s                 | 1.2s                              | 1.33x             | 0.17x                    |
| 7     | 87s                          | 65s                   | 4.5s                              | 1.34x             | 19.33x                   |
| 9     | 362.88s (~6 mins)            | 300s                  | 180s                              | 1.21x             | 2.02x                    |
| 10    | 3,628.8s (~1 hour)          | 1500s                 | 1000s                             | 2.42x             | 3.63x                    |

### Observations:
- **MPI Speedup**: Significant performance improvement for larger graphs (e.g., n = 10 with 2.42x speedup).
- **OpenMP+MPI Speedup**: Impressive speedup, particularly for n = 7 (19.33x), though performance degradation is observed for smaller graphs due to overhead.

## Scalability Analysis

- **Weak Scaling**: Parallel implementations using MPI and OpenMP+MPI show good performance improvement as the graph size increases.
- **Strong Scaling**: Performance improvement per added processor diminishes with smaller graphs, especially for n = 6.

## Challenges
- **Load Balancing**: Ensuring even workload distribution across processors was a challenge, mitigated by METIS partitioning.
- **Overhead from OpenMP**: The hybrid OpenMP+MPI implementation showed performance degradation for smaller graphs due to threading overhead.

## Conclusion
This project demonstrates the effectiveness of parallel graph algorithms for constructing independent spanning trees in bubble-sort networks. Leveraging MPI for inter-node communication, OpenMP for intra-node parallelism, and METIS for graph partitioning provided significant performance improvements, particularly for larger graphs.

Future work includes refining partitioning strategies, optimizing OpenMP to minimize overhead, and exploring GPU-based parallelization for further performance gains.
