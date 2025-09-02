# CA-FinalProject-RNA-Quant-RISCV
Computer Architecture final project - Accelerate RNA Sequence Quantification on RISC-V Distributed Cluster on Near-DRAM Platform

### Project Theme
- Perform RNA sequence quantification by matching query sequences (R) to a set of transcripts (T).
- Build a k-mer index of transcripts and use it to find the transcript that each query sequence is most similar to.
- Distribute the workload across multiple logic units to simulate parallel processing on a near-DRAM RISC-V platform.
- **Final output:** For each transcript T, count how many queries R match it best (e.g., `T0: 2, T1: 4, T3: 0, T4: ...`).
- **Objective:** Minimize execution time and maximize performance by efficiently distributing and processing sequences in parallel.

## Folder Structure
- `project_spec/` : Original project instructions and PPT provided by the TA
- `report/` : Final report including analysis and code explanations
- `src/` : Source code
- `README.md` : Project description

## Implementation Details
1. **Workload Distribution (`distribute.cpp`)**
   - Simulate distributing transcripts (T) and query sequences (R) across 256 logic units
   - Split the sequences evenly among units to simulate parallel processing

2. **Indexing (`index.c` & `hash-map.h`)**
   - Extract k-mers of a specified length from transcripts (T)
   - Build a hashmap that maps each k-mer to the transcripts where it appears
   - **Hashmap implementation (`hash-map.h`):** Contains the data structures and functions for creating and querying the k-mer hashmap efficiently

3. **Quantification (`quantify.c` & `hash-map.h`)**
   - Extract k-mers from query sequences (R) of the same length as in the index
   - Use functions in `hash-map.h` to query the hashmap and find which transcript each query matches best
   - **Final output:** For each transcript T, count how many queries R match it best  
     - Example: `T0: 2, T1: 4, T3: 0, T4: ...`

4. **Simulation (`decoder.isa` & gem5)**
   - Run the system in a **gem5 simulated RISC-V environment**
   - **Custom ISA modifications (`decoder.isa`):** Implemented custom instructions for the Accelerator and Vector Coprocessor to speed up hash table operations and parallel k-mer processing
