# **BinarizeDynamics.jl Specification**

Version: 1.0-draft  
Last Updated: 2025-12-30  
Status: Draft  
Language: Julia 1.10+

## **1\. Project Overview & Environment Context**

### **1.1 Package Identity**

|

| Field | Value |  
| Package Name | BinarizeDynamics.jl |  
| One-line Description | High-performance analysis of structural dynamics in sequence data via pairwise binarization and bitwise operations. |  
| Language | Julia 1.10+ |  
| License | MIT |  
| Author | Project Maintainer |  
| Repository | TBD |

### **1.2 Environment Context**

| Requirement | Specification |  
| Target Language Version | Julia 1.10+ (LTS or latest stable) |  
| Dependency Manager | Pkg (Julia built-in) |  
| Configuration Files | Project.toml, Manifest.toml |  
| OS Compatibility | Linux, macOS, Windows |  
| Virtual Environment | Project.toml environment activation |

## **2\. Problem Statement**

### **2.1 Current State (As-Is)**

**PROBLEM**: Analyzing structural linkages and dynamic changes between groups of aligned sequences (DNA, logs, etc.) is computationally prohibitive using standard character comparisons, leading to memory bottlenecks (O(N²) char storage) and slow execution speeds.

**CURRENT\_SOLUTIONS**:

* **Standard MSA tools (BioPython/BioJulia)**: Focus on static consensus or simple conservation, missing "correlated mutation" dynamics.  
* **BinarizeComp.jl**: Calculates mismatch rates within positions but discards the specific "pair" information (i, j) required to track *who* changed with *whom* across different positions.

**PAIN\_POINTS**:

1. **Memory Explosion**: Storing full character matrices for all pairwise comparisons.  
2. **Missing Linkage**: Inability to trace if the specific pair (i, j) that differs at Position A is the *same* pair differing at Position B.

### **2.2 Desired State (To-Be)**

**SOLUTION**: A specialized Julia package that abstracts sequence relationships into binary (0/1) signals using BitArray storage, enabling high-speed 2D structural analysis and 2.5D differential dynamics.

**KEY\_IMPROVEMENTS**:

1. **Memory Efficiency**: 8x reduction using BitArray vs. Byte storage.  
2. **Speed**: 10\~100x speedup via bitwise operations (XOR, AND, Popcount) and threading.  
3. **Deep Structural Insight**: Reveals "coupling" (positive correlation) and "decoupling" (negative correlation) between positions.

### **2.3 Target Users**

| User Type | Technical Level | Primary Use Case |  
| Bioinformatician | Advanced | Protein/DNA co-evolution and structural constraint analysis. |  
| Data Scientist | Intermediate | Analyzing coupled failure modes in system logs or financial tickers. |  
| System Engineer | Advanced | Tracing error propagation paths in distributed systems. |

### **2.4 Competitive Analysis**

| Package | Language | Limitation vs This Package |  
| BitMatrix.jl | Julia | Name conflict; generic bit tools only, no domain logic. |  
| BinarizeComp.jl | Julia | Aggregates stats per position; loses pairwise identity linkage. |  
| Evoformer (AF2) | Python/JAX | Embedded in massive DL framework; not usable as standalone analysis tool. |  
**UNIQUE\_VALUE\_PROPOSITION**: Standalone, lightweight tool that democratizes "Evoformer-like" structural pair analysis using efficient bitwise mathematics without heavy deep learning dependencies.

## **3\. Architecture & Modularization Strategy**

### **3.1 Root Directory Structure**

BinarizeDynamics.jl/  
├── src/  
│   ├── BinarizeDynamics.jl       \# Main module, re-exports  
│   ├── types.jl                  \# Structs: PairMapper, BinarizedData  
│   ├── utils.jl                  \# Bitwise helpers, validation  
│   ├── binarization.jl           \# Core: Ingest \-\> 3D BitArray  
│   ├── analysis.jl               \# Core: Phi Coeff, APC, Differential  
│   └── visualization.jl          \# Plotting recipes (Makie extension)  
├── test/  
│   ├── runtests.jl  
│   ├── test\_types.jl  
│   ├── test\_binarization.jl  
│   ├── test\_analysis.jl  
│   └── fixtures/                 \# Mock data generators  
├── docs/  
│   ├── make.jl  
│   └── src/  
├── Project.toml  
├── Manifest.toml  
├── README.md  
└── .gitignore

### **3.2 Modularization Principles**

| Principle | Implementation |  
| Separation of Concerns | binarization.jl handles data transformation; analysis.jl handles statistics; visualization.jl handles plotting. |  
| High Cohesion | All bitwise math logic resides in utils.jl or localized core functions to allow SIMD optimization. |  
| Scalability | Designed to support a SparseArrays extension in the future for massive datasets. |

### **3.3 Module Responsibilities**

| Module | Responsibility | Exports | Dependencies |  
| types.jl | Define custom data structures | PairMapper, BinarizedData, InteractionMatrix | None |  
| utils.jl | Helper functions for indexing and validation | validate\_sequences, index\_to\_pair | types |  
| binarization.jl | Convert raw strings to 3D binary structures | binarize | types, utils |  
| analysis.jl | Statistical analysis (Phi, APC, Diff) | analyze\_structure, diff\_dynamics | types, utils |  
| visualization.jl | Plotting recipes | plot\_interaction, plot\_differential | types, Makie |

### **3.4 Implementation Order**

1\. types.jl          (Define PairMapper, BinarizedData)  
2\. utils.jl          (Implement index mapping logic binom(n,2))  
3\. binarization.jl   (Implement core 3D binarization loop)  
4\. analysis.jl       (Implement Phi coefficient and APC)  
5\. BinarizeDynamics.jl (Main entry point)  
6\. visualization.jl  (Optional/Extension)

## **4\. Functional Requirements & Self-Verification**

### **4.1 MUST Have Requirements**

| ID | Function | Description | Input | Output | Acceptance Criteria |  
| FR-001 | binarize | Convert N sequences of length L into pairwise binary differences. | Vector{String} (length N) | BinarizedData (L × Pairs × 1\) | BitArray dimensions correct; pair mapping accurate. |  
| FR-002 | analyze\_structure | Calculate structural linkage (Phi coefficient) \+ APC correction. | BinarizedData | InteractionMatrix (L × L) | Returns matrix with values  
$$-1, 1$$  
; diagonal is 1.0. |  
| FR-003 | diff\_dynamics | Compute differential dynamics between two interaction matrices. | InteractionMatrix (A), InteractionMatrix (B) | InteractionMatrix (Diff) | Matrix\_B \- Matrix\_A; background noise reduced. |

### **4.2 Detailed Requirement Specifications**

**FR-001: binarize**

INPUT SPECIFICATION:  
  sequences:  
    type: Vector{String}  
    constraints:  
      \- Length \> 1  
      \- All strings must have equal length L  
    optional: false  
  nthreads:  
    type: Int  
    default: Threads.nthreads()

OUTPUT SPECIFICATION:  
  type: BinarizedData  
  structure:  
    data: BitArray{3} (dims: L positions × P pairs × 1\)  
    mapper: PairMapper (handles (i,j) \<-\> index logic)  
  guarantees:  
    \- P \== binomial(N, 2\)  
    \- data\[pos, idx\] \== 1 if seqs\[i\]\[pos\] \!= seqs\[j\]\[pos\]

**FR-002: analyze\_structure**

INPUT SPECIFICATION:  
  bin\_data:  
    type: BinarizedData  
  method:  
    type: Symbol  
    default: :phi\_apc  
    valid\_values: \[:phi, :phi\_apc\]

OUTPUT SPECIFICATION:  
  type: InteractionMatrix  
  structure:  
    matrix: Matrix{Float64} (L × L)  
    positions: Vector{Int}  
  guarantees:  
    \- Symmetric matrix  
    \- Values represent correlation (Phi coefficient)

### **4.3 Self-Verification Protocol**

**AI SELF-AUDIT INSTRUCTION**:

After generating code, verify:

IMPLEMENTATION\_VERIFICATION:  
  FR-001 (binarize):  
    \- \[ \] Input sequence alignment verified (all lengths equal)  
    \- \[ \] BitArray dimensions match (L, n\*(n-1)/2, 1\) or (L, Pairs) depending on optimized layout  
    \- \[ \] PairMapper correctly converts (i,j) tuples to linear indices  
    
  FR-002 (analyze\_structure):  
    \- \[ \] Phi calculation handles divide-by-zero (e.g. constant columns)  
    \- \[ \] APC correction formula applied correctly: Raw \- (Mean\_i \* Mean\_j / Mean\_all)  
    \- \[ \] Output matrix is L x L symmetric

## **5\. Technical Specifications**

### **5.1 Runtime Environment**

| Requirement | Specification |  
| Language | Julia 1.10+ |  
| OS | Linux, macOS, Windows |  
| Minimum RAM | 16GB (recommended for N\>1000, L\>1000) |  
| Parallelism | Multi-threading enabled (Threads.@threads) |

### **5.2 Environment Setup**

\# Julia REPL  
using Pkg  
Pkg.generate("BinarizeDynamics")  
Pkg.activate("BinarizeDynamics")  
Pkg.add("Statistics")  
Pkg.add("LinearAlgebra")  
Pkg.add("Makie") \# For visualization extension

### **5.3 Dependencies**

#### **Core Dependencies**

| Package | Version | Purpose | Justification |  
| Statistics | StdLib | Mean, Var calc | Core stats |  
| LinearAlgebra | StdLib | Matrix ops | Matrix subtractions, diagonals |  
| LoopVectorization | \>=0.12 | SIMD optimization | Accelerate bit counting/stats |

#### **Development Dependencies**

| Package | Version | Purpose |  
| Test | StdLib | Unit testing |  
| BenchmarkTools | \>=1.3 | Performance profiling |

### **5.4 Coding Standards**

| Element | Convention | Example |  
| Functions | snake\_case | analyze\_structure (per Blueprint & Spec) |  
| Types | PascalCase | BinarizedData |  
| Variable | snake\_case | pair\_index |  
| Formatter | JuliaFormatter | style \= "blue" |

## **6\. API Design**

### **6.1 Type Definitions**

#### **PairMapper**

TYPE: PairMapper  
CATEGORY: struct  
DESCRIPTION: Maps string indices (i, j) to linear pair index k.

FIELDS:  
  n\_sequences:  
    type: Int  
    description: Number of original sequences (N)  
  total\_pairs:  
    type: Int  
    description: binomial(N, 2\)

CONSTRUCTORS:  
  \- PairMapper(n::Int): Calculates total\_pairs automatically.

#### **BinarizedData**

TYPE: BinarizedData  
CATEGORY: struct  
DESCRIPTION: Container for 3D binary difference data.

FIELDS:  
  data:  
    type: BitMatrix  
    description: Dimensions (Pairs × Positions). Note: Transposed for memory locality if needed.  
  mapper:  
    type: PairMapper  
    description: Associated mapping logic.  
  seq\_length:  
    type: Int  
    description: L

INVARIANTS:  
  \- size(data, 1\) \== mapper.total\_pairs  
  \- size(data, 2\) \== seq\_length

#### **InteractionMatrix**

TYPE: InteractionMatrix  
CATEGORY: struct  
DESCRIPTION: Result of 2D structural analysis.

FIELDS:  
  matrix:  
    type: Matrix{Float64}  
    description: L × L correlation matrix.  
  method:  
    type: Symbol  
    description: Method used (:phi, :phi\_apc).

### **6.2 Function Definitions**

#### **binarize**

FUNCTION: binarize  
SIGNATURE: binarize(sequences::Vector{String}; nthreads::Int=Threads.nthreads()) \-\> BinarizedData  
IMPLEMENTS: FR-001

PURPOSE: Compare all pairs of sequences and generate binary difference map.

PARAMETERS:  
  sequences:  
    type: Vector{String}  
    description: Aligned sequences.

RETURNS:  
  type: BinarizedData  
  description: The compressed binary representation.

COMPLEXITY: O(L \* N^2) time, O(L \* N^2) bits space.

#### **analyze\_structure**

FUNCTION: analyze\_structure  
SIGNATURE: analyze\_structure(data::BinarizedData; method::Symbol=:phi\_apc) \-\> InteractionMatrix  
IMPLEMENTS: FR-002

PURPOSE: Compute structural correlations between positions.

PARAMETERS:  
  data:  
    type: BinarizedData  
  method:  
    type: Symbol  
    valid\_values: :phi, :phi\_apc (Average Product Correction)

RETURNS:  
  type: InteractionMatrix  
  guarantees: values in range \[-1.0, 1.0\] usually (APC can shift this).

#### **diff\_dynamics**

FUNCTION: diff\_dynamics  
SIGNATURE: diff\_dynamics(mat\_b::InteractionMatrix, mat\_a::InteractionMatrix) \-\> InteractionMatrix  
IMPLEMENTS: FR-003

PURPOSE: Calculate the differential dynamics (Matrix B \- Matrix A).

ERRORS:  
  DimensionMismatch: If matrices have different sizes.

## **7\. Workflow & Data Pipeline**

### **7.1 User Workflow: Structural Analysis**

WORKFLOW\_ID: WF-001 IMPLEMENTS: FR-001, FR-002, FR-004

#### **Steps**

STEP 1: Load & Ingest  
  INPUT: raw\_sequences (Vector{String})  
  PROCESS: binarize(raw\_sequences)  
  OUTPUT: bin\_data (BinarizedData)  
  NEXT: STEP 2

STEP 2: Structural Analysis  
  INPUT: bin\_data  
  PROCESS: analyze\_structure(bin\_data, method=:phi\_apc)  
  OUTPUT: interaction\_matrix (InteractionMatrix)  
  NEXT: STEP 3

STEP 3: Visualization  
  INPUT: interaction\_matrix  
  PROCESS: plot\_interaction(interaction\_matrix)  
  OUTPUT: Plot object (Makie)

### **7.2 Data Flow**

Raw Strings (Vector{String})  
│ Size: N \* L bytes  
│  
▼ ──── binarize() ────  
│  
BinarizedData (BitMatrix)  
│ Size: (N\*(N-1)/2 \* L) / 8 bytes  
│ Reduction: \~8x vs storing differences as bytes  
│  
▼ ──── analyze\_structure() ────  
│  
InteractionMatrix (Matrix{Float64})  
│ Size: L \* L \* 8 bytes  
│ Abstracted structure

## **8\. Data Specifications**

### **8.1 Input Formats**

#### **Aligned Sequences (Memory)**

TYPE: Vector{String}  
CONSTRAINT:  
  \- All strings equal length.  
  \- No empty vector.

### **8.2 Output Formats**

#### **Interaction Heatmap**

| Field | Type | Description |  
| x-axis | Int | Position Index 1 |  
| y-axis | Int | Position Index 2 |  
| color | Float | Phi/APC value (Intensity of linkage) |

## **9\. Testing Strategy**

### **9.1 Test Categories**

| Category | Scope | Tool |  
| Unit | PairMapper index logic, Phi calc math | Test |  
| Integration | End-to-end binarize \-\> analyze | Test |  
| Performance | Benchmarking large N | BenchmarkTools |

### **9.2 Test Cases**

#### **FR-001: binarize**

| TC-ID | Description | Input | Expected |  
| TC-01 | Minimal Case | \["A", "B"\] (L=1) | 1 pair, value=1 (diff) |  
| TC-02 | Identical | \["A", "A"\] | 1 pair, value=0 (same) |  
| TC-03 | Length Mismatch | \["A", "AA"\] | ArgumentError |

#### **FR-002: analyze\_structure (Simulation)**

| TC-ID | Description | Setup | Expected |  
| TC-SIM-01 | Perfect Coupling | Pos 1 & 2 flip together | High correlation (\>0.9) |  
| TC-SIM-02 | Random Noise | Random sequences | Low correlation (\~0.0) |

## **10\. Documentation Plan**

### **10.1 Document Types**

| Document | Audience | Location |  
| README | Users | README.md (Quickstart, Install) |  
| API Ref | Devs | docs/src/api.md |  
| Tutorial | Users | docs/src/tutorial.md (The Simulation Workflow) |

### **10.2 Docstring Standard**

Follow Julia Docstring standard:

"""  
    binarize(sequences::Vector{String}) \-\> BinarizedData

Convert sequences into pairwise binary difference maps.

\# Arguments  
\- \`sequences\`: Aligned sequence strings.

\# Returns  
\- \`BinarizedData\`: Struct containing the BitMatrix.  
"""  
function binarize...

## **11\. Development Timeline**

### **11.1 Phases**

| Phase | Name | Focus | Deliverables |  
| 1 | Skeleton | Proj Structure, Types | types.jl, utils.jl |  
| 2 | Core MVP | Binarization Logic | binarization.jl, tests |  
| 3 | Analysis | Stats & Math | analysis.jl (Phi/APC) |  
| 4 | Visuals | Plots & Docs | visualization.jl, Documentation |

## **12\. Version Control & CI/CD**

### **12.1 Git Setup**

git init  
git add .  
git commit \-m "Initial commit: BinarizeDynamics skeleton"

### **12.2 CI Workflow (GitHub Actions)**

* **Trigger**: Push to main, PRs.  
* **Job**: Run julia \--project \-e 'using Pkg; Pkg.test()'  
* **Matrix**: Julia 1.10, 1.9, OS (Ubuntu, Windows).

## **Specification Checklist**

* $$x$$  
  All 12 sections present  
* $$x$$  
  Package Name: BinarizeDynamics.jl  
* $$x$$  
  Language: Julia  
* $$x$$  
  Function names: snake\_case (per blueprint/spec rules)  
* $$x$$  
  Core Logic: Pairwise Binarization \+ Phi/APC  
* $$x$$  
  Self-Verification Protocols Included