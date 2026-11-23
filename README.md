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
- [Containers](#containers)
- [Quality Assurance](#quality-assurance)
- [Acknowledgements](#acknowledgements)
- [References](#references)
- [Citation](#citation)
- [License](#license)

## Key Files

| File                              | Notes                                                                                      |
| :-------------------------------- | :----------------------------------------------------------------------------------------- |
| `src/mpi-greedy-set-cover.cpp`    | C++ program.                                                                               |
| `data/coordinates.json`<br>&nbsp; | Test data (3,000 targets, 124 KB).<br>&nbsp;&nbsp;&nbsp;Galactic latitudes and longitudes. |
| `apptainer.def`                   | Optional. Apptainer definition file.                                                       |

## Software Requirements

| Software     | Notes                                                                                                                      |
| :----------- | :------------------------------------------------------------------------------------------------------------------------- |
| C++          | [Details here](https://isocpp.org/get-started). Free and proprietary compilers available.                                  |
| JSON for C++ | <!-- textlint-disable terminology -->[nlohmann/json](https://json.nlohmann.me/)<!-- textlint-enable terminology -->. Free. |
| MPI          | Multiple implementations available. Free.                                                                                  |
| Apptainer    | Optional. [Available here](https://apptainer.org/). Free.                                                                  |

## Getting Started

### Compilation

Program source code requires compilation on the target system with a C++ / MPI
compiler. E.g., using the mpicxx wrapper with G++, the GNU C++ compiler,:

    mpicxx src/mpi-greedy-set-cover.cpp -o mpi-greedy-set-cover.exe -O3

#### G++ Compiler Flags

| Compiler Flag | Notes                              |
| :------------ | :--------------------------------- |
| `-o`          | Name of the output file.           |
| `-O3`         | Optional. Aggressive optimisation. |

### Execution

The compiled C++ program may be executed from a command line. A JSON data file
to process should be specified as an input argument.

Execution syntax for some example command lines, the MSYS2 UCRT64 and Ubuntu
shells, is given below. Note, an example 3,000-target JSON file has been
supplied in this repository.

    mpiexec -n 4 mpi-greedy-set-cover.exe data/coordinates.json

#### Command Line Options

| Option         | Notes                                           |
| :------------- | :---------------------------------------------- |
| `-n <num>`     | Number of processes to allocate to the program. |
| `<input_file>` | A JSON data file to process.                    |

#### Open MPI in WSL

Open MPI in the Windows Subsystem for Linux (WSL) doesn't always recognise all
of the processes a Linux distribution exposes. Invocation with an option
directing use of hardware threads as valid processing slots may help:

    mpiexec --use-hwthread-cpus -n 4 mpi-greedy-set-cover.exe data/coordinates.json

## Containers

An optional Apptainer container definition file has been provided. The build
process compiles the set cover program to `/opt/app/mpi-greedy-set-cover.exe`
in an Ubuntu container. Example data is built into
`/opt/app/data/coordinates.json`.

Assuming execution from the repository root, example Ubuntu container build
syntax is:

    sudo apptainer build mpi-greedy-set-cover.sif apptainer.def

If Apptainer will be invoked in rootless mode, Linux may need the `squashfuse`
package to mount components. In Ubuntu, the package may be installed via:

    sudo apt install squashfuse

Example containerised invocation syntax is given below, for four processes.

    mpiexec -n 4 apptainer exec mpi-greedy-set-cover.sif /opt/app/mpi-greedy-set-cover.exe /opt/app/data/coordinates.json

Example containerised invocation, addressing Open MPI / WSL process recognition
issues, is given below.

    mpiexec --use-hwthread-cpus -n 4 apptainer exec mpi-greedy-set-cover.sif /opt/app/mpi-greedy-set-cover.exe /opt/app/data/coordinates.json

## Quality Assurance

The repository code has been tested across 1, 2, 4 and 8 processes, in the
following environments.

<details>
<summary>Ubuntu 26.04, G++ 15.2.0</summary>

<br>

| Type     | Component                                                                              | Version            |
| :------- | :------------------------------------------------------------------------------------- | :----------------- |
| Platform | Operating system                                                                       | Ubuntu 26.04       |
| Software | GNU C++                                                                                | 15.2.0             |
| &quot;   | Open MPI                                                                               | 5.0.10             |
| &quot;   | Apptainer                                                                              | 1.4.5              |
| Packages | <!-- textlint-disable terminology -->nlohmann/json<!-- textlint-enable terminology --> | 3.12.0             |
| Data     | `coordinates.json`                                                                     | Repository dataset |

</details>

<!-- jscpd:ignore-start -->

<details>
<summary>Windows 11, G++ 16.2.0</summary>

<br>

| Type     | Component                                                                              | Version                                               |
| :------- | :------------------------------------------------------------------------------------- | :---------------------------------------------------- |
| Platform | Operating system                                                                       | Windows 11, 26H2 (OS Build 26300.9818)                |
| Software | GNU C++                                                                                | 16.2.0, run from MSYS2 version 20260927, UCRT64 shell |
| &quot;   | Microsoft MPI                                                                          | 10.1.12498.52                                         |
| Packages | <!-- textlint-disable terminology -->nlohmann/json<!-- textlint-enable terminology --> | 3.11.3                                                |
| Data     | `coordinates.json`                                                                     | Repository dataset                                    |

</details>

<!-- jscpd:ignore-end -->

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

3. _Microsoft Copilot_. (2026). Microsoft. Accessed: 29 October to 07 November, 2026.
   [AI model output, Microsoft Copilot application (cloud‑based service)]
   Available:
   [Download the Microsoft Copilot app](https://www.microsoft.com/en-us/copilot/download-copilot-app).

## Citation

Citation details are available by clicking **Cite this repository** in the
GitHub sidebar.

## License

This repository is licensed under the BSD-3-Clause license. Details are
available in the [LICENSE](./LICENSE) file.
