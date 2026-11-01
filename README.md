# Greedy Set Cover Field Selection for Multi-object Spectroscopy in C++ MPI

[![super-linter](../../actions/workflows/super-linter.yml/badge.svg)](../../actions/workflows/super-linter.yml)

This repository contains code and supporting material associated with the
article "Greedy Set Cover Field Selection for Multi-object Spectroscopy in C++
MPI" [[1](#references)]. The article explored the problem of efficiently
observing large populations of astronomical targets that cannot be covered by a
single telescope pointing. Using a greedy set-cover approach, the software
repeatedly selects the field of view that captures the greatest number of
previously unobserved targets, maximising the scientific return obtainable from
a limited number of telescope pointings.

---

<figure>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="assets/greedy-set-cover-object-pairs-dm.svg">
    <img src="assets/greedy-set-cover-object-pairs-lm.svg" loading="lazy" alt="Scattered points within three circular regions; two overlap and one is separate, with selected points in red." width="100%">
  </picture>
  <br>
  <figcaption>Figure 1. Example set covers. For cover radius <i>r</i> and object distance <i>d</i>, example covers for red object pairs at the <i>d</i> = 2<i>r</i> and <i>d</i> < 2<i>r</i> cases are shown.</figcaption>
</figure>

---

## Table of Contents

- [Key Files](#key-files)
- [Software Requirements](#software-requirements)
- [Getting Started](#getting-started)
- [Quality Assurance](#quality-assurance)
- [Acknowledgements](#acknowledgements)
- [References](#references)
- [Citation](#citation)
- [License](#license)

## Key Files

| File                           | Notes        |
| :----------------------------- | :----------- |
| `src/mpi-greedy-set-cover.cpp` | C++ program. |

## Software Requirements

| Software | Notes                                                                                     |
| :------- | :---------------------------------------------------------------------------------------- |
| C++      | [Details here](https://isocpp.org/get-started). Free and proprietary compilers available. |
| MPI      | Free. Multiple implementations available.                                                 |

## Getting Started

### Compilation

Program source code requires compilation on the target system with a C++ / MPI
compiler. E.g., using the mpicxx wrapper with G++, the GNU C++ compiler,:

    mpicxx src/mpi-greedy-set-cover.cpp -o mpi-greedy-set-cover.exe

#### G++ Compiler Flags

| Compiler Flag | Notes                    |
| :------------ | :----------------------- |
| `-o`          | Name of the output file. |

### Execution

The compiled C++ program should be executed from a command line. Execution
syntax for an example command line, the MSYS2 UCRT64 shell, is given below.
Execution over four processes was specified.

    mpiexec -n 4 mpi-greedy-set-cover.exe

## Quality Assurance

The repository code has been tested in the following environment.

<details>
<summary>Windows 11, G++ 16.2.0</summary>

<br>

| Type     | Component        | Version                                               |
| :------- | :--------------- | :---------------------------------------------------- |
| Platform | Operating system | Windows 11, 26H2 (OS Build 26300.9818)                |
| Software | GNU C++          | 16.2.0, run from MSYS2 version 20260927, UCRT64 shell |
| &quot;   | Microsoft MPI    | 10.1.12498.52                                         |

</details>

## Acknowledgements

This work was supported by a Macquarie University Research Excellence
Scholarship. Duck.ai [[2](#references)] and Microsoft Copilot [[3](#references)]
were used as software engineering and repository management aids.

## References

1. T. N. Stenborg, "Greedy Set Cover Field Selection for Multi-object
   Spectroscopy in C++ MPI", in _Astron. Data Anal. Softw. Syst. XXIV_, in
   Astronomical Society of the Pacific Conference Series, vol. 495, A. R.
   Taylor and E. Rosolowsky, Eds., 2015, pp. 269&ndash;272.\
   [View PDF](https://aspbooks.org/publications/495/269.pdf)
   &nbsp; [View at publisher](https://aspbooks.org/custom/publications/paper/495-0269.html)
   &nbsp; [SciX](https://scixplorer.org/abs/2015ASPC..495..269S/abstract)

2. _Duck.ai_. (Large language model, GPT-6 Luna). DuckDuckGo.
   [Online]. Available: [duck.ai](https://duck.ai/).

3. _Microsoft Copilot_. (2026). Microsoft. Accessed: 29 October to 01 November, 2026.
   [AI model output, Microsoft Copilot application (cloud‑based service)]
   Available:
   [Download the Microsoft Copilot app](https://www.microsoft.com/en-us/copilot/download-copilot-app).

## Citation

Citation details are available by clicking **Cite this repository** in the
GitHub sidebar.

## License

This repository is licensed under the BSD-3-Clause license. Details are
available in the [LICENSE](./LICENSE) file.
