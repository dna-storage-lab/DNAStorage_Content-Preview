
<h1 id="TJU-RuiQin">Content-Derived Molecular Access for On-Demand Preview in DNA Data Storage</h1>

## Table of Contents

- [About the Project](#about-the-project)
- [Project Structure](#project-structure)
- [Source Data](#source-data)
- [Pipeline](#pipeline)
- [Dependencies](#dependencies)
- [License](#license)

## About the Project

Synthetic DNA is a highly promising medium for high-density data storage. However, most current archival DNA storage systems still rely on whole-library sequencing and computational reconstruction, while rapid content-based partial readout remains challenging. This repository provides the computational reconstruction pipeline for processing nanopore sequencing data generated after content-derived PCR amplification. The molecular workflow first enriches selected target regions using content-derived primer pairs. The software then processes the resulting basecalled FASTQ reads together with the corresponding ideal companion sequences to recover the previewed content.

The implemented pipeline performs repeated test sampling, nucleotide demapping, companion-sequence-guided target-read alignment, forward–backward probabilistic inference, multi-read probability merging, and final data reconstruction. It supports the eight target-content preview experiments included in this project.

This repository focuses on post-sequencing read processing and data reconstruction.

---

## Project Structure


```text
/
├── companion_sequence/                 # Ideal companion sequences
│   ├── companion_sequence_poem1.txt
│   ├── companion_sequence_poem2.txt
│   ├── companion_sequence_poem3.txt
│   ├── companion_sequence_poem4.txt
│   ├── companion_sequence_poem5.txt
│   ├── companion_sequence_poem6.txt
│   ├── companion_sequence_poem7.txt
│   └── companion_sequence_poem8.txt
├── sequencing_data/
│   ├── poem1.fastq
│   ├── poem2.fastq
│   ├── poem3.fastq
│   ├── poem4&5.fastq
│   ├── poem6.fastq
│   ├── poem7.fastq
│   └── poem8.fastq
├── src/                                    # C/C++ sources and Bash scripts
│   ├── read_sampling.cpp                   # Implementation of repeated FASTQ sampling
│   ├── run_read_sampling.sh
│   ├── edlib.h
│   ├── edlib.cpp
│   ├── companion_sequence_alignment.cpp    # Implementation of demapping and alignment
│   ├── run_companion_sequence_alignment.sh
│   ├── HalfMaker.h
│   ├── fb_indel_correction.cpp             # Implementation of forward–backward inference
│   ├── fb_probability_extraction.cpp       # Implementation of probability extraction
│   ├── run_fb_indel_correction.sh
│   ├── probability_merging.cpp             # Implementation of probability merging
│   ├── preview_data_recovery.cpp           # Implementation of content reconstruction
│   └── run_preview_data_recovery.sh
└──README.md
```

## Source Data

Raw nanopore sequencing data are available in the NCBI Sequence Read Archive under accession number PRJNA1505625. The provided sequencing datasets were generated using a PromethION 2 Solo nanopore sequencer and basecalled with the super-accuracy (SUP) model.

## Pipeline

### 1. Repeated FASTQ read sampling

Command:

```bash
bash src/run_read_sampling.sh
```

- **Input:** Basecalled FASTQ files in `Content_preview/sequencing_data/`.
- **Processing:** The program reads complete four-line FASTQ records and creates consecutive, non-overlapping repeated-test subsets. 
- **Output:** `result/sampled_fastq/<dataset>/coverage_<N>/replicate_<RRR>.fastq`.

The word `coverage` is retained in directory names for compatibility with the original analysis. Here, `N` is the number of raw FASTQ reads assigned to one replicate.

### 2. Demapping and companion-sequence alignment

Command:

```bash
bash src/run_companion_sequence_alignment.sh
```

- **Input:** Sampled FASTQ files from Step 1 and the corresponding ideal sequence in `Content_preview/companion_sequence/`.
- **Processing:** Each nucleotide is demapped according to `{A→00, T→01, G→10, C→11}`. The upper bit forms the noisy information sequence and the lower bit forms the noisy companion sequence. Edlib aligns the noisy companion sequence to its ideal companion sequence and retains target-matched reads with positional information.
- **Output:** `result/aligned_reads/<target>/coverage_<N>/replicate_<RRR>_match_reads_lower.txt`.

### 3. Forward–backward probabilistic inference

Command:

```bash
bash src/run_fb_indel_correction.sh
```

- **Input:** Target-matched reads from Step 2.
- **Processing:** A forward–backward model accounts for insertion and deletion events and estimates the probability of binary states 0 and 1 at each information-bit position.
- **Output:** Corrected-read, `prob0`, and `prob1` files under `result/fb_results/`, preserving the target, coverage, and replicate hierarchy.

### 4. Probability merging and preview-data reconstruction

Command:

```bash
bash src/run_preview_data_recovery.sh
```

- **Input:** Position-wise probability files from Step 3.
- **Processing:** Log probabilities are accumulated across the aligned reads in each Step 1 replicate. Maximum-log-likelihood decisions recover the information bits, which are then converted into the preview content.
- **Output:** Voting results and reconstructed preview files under `result/recovered/`.


## Dependencies

- Linux with Bash
- `g++` with C++11 and C++17 support
- POSIX `find` and `sort`
- Edlib source and API header (included as `src/edlib.cpp` and `src/edlib.h`)

## License

This project is licensed under the MIT License. Edlib remains subject to its upstream license.