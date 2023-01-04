---
title: 'excitingtools: An exciting Workflow Tool'
tags:
 - Python
 - Physics
 - DFT
 - Workflow
authors:
  - name: Alexander Buccheri
    orcid: 0000-0001-5983-8631
    affiliation: 1, 2
  - name: Fabian Peschel
    orcid: 0000-0003-0619-6713
    affiliation: 1
  - name: Benedikt Maurer
    orcid: 0000-0001-9152-7390
    affiliation: 1
  - name: Mara Voiculescu
    orcid: 0000-0003-4393-8528
    affiliation: 1
  - name: Daniel T. Speckhard
    orcid: 0000-0002-9849-0022
    affiliation: 1
  - name: Hannah Kleine
    orcid: 0000-0003-2251-8719
    affiliation: 1
  - name: Elisa Stephan
    orcid: 0000-0002-6359-9044
    affiliation: 1
  - name: Martin Kuban
    orcid: 0000-0002-1619-2460
    affiliation: 1
  - name: Claudia Draxl
    orcid: 0000-0003-3523-6657
    affiliation: 1
affiliations:
- name: Humboldt-Universität zu Berlin, Berlin, Germany
  index: 1
- name: Max Planck Institute for the Structure and Dynamics of Matter, Hamburg, Germany
  index: 2
date: December 2022
bibliography: paper.bib
---

# Introduction

<span style="font-family:american typewriter; font-size:1em;">**exciting**</span> [@gulans:2014] is a full-potential, 
all-electron, density-functional theory (DFT) code, which uses linearized augmented plane waves plus local orbitals 
(LAPWs and LOs, respectively) for valence and semi-core electrons, and explicitly treats core electrons via the 
Dirac equation. This basis set provides a systematic path for reaching the complete-basis-set limit, relying only on 
well-controlled numerical approximations. The high precision of the LAPW basis set makes <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> 
well-suited for the generation of benchmark-quality data, and to serve as a reference for other DFT implementations, 
especially those relying on pseudopotentials.

It has been shown that <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> is able to 
achieve microhartree precision for total energies in DFT calculations [@gulans:2018]. Also in G~0~W~0~ calculations, the complete-basis-set limit limit was attained [@nabok:2016]. Its high precision has also been 
unequivocally demonstrated in the so-called $\Delta$ test [@lejaeghere:2016], which compared the relative precision of several 
DFT codes for a benchmark set of 71 elemental crystals.

# Statement of Need

Materials databases such as NOMAD [@draxl:2019], AFLOW [@curtarolo:2012], and the Materials Project [@jain:2013] 
host millions of DFT results. The majority of data have been computed with pseudopotential-based DFT codes, and is thus lacking validation with more precise methods. There is a strong need to provide databases with benchmark-quality results, 
which serve to give an indication of the precision one can achieve in a given material property, with a specific method 
and settings. For scientists and engineers who wish to compute specific properties to some required precision, having an 
indication of the optimal settings and suitable DFT approximations is extremely valuable. Beyond the ground state, 
materials databases for excited state calculations are, in general, strongly lacking. The generation of large amounts of excited state data will require both reliable ground state calculations as inputs and analogous benchmark-quality calculations. Moreover, machine learning models that predict material properties would greatly benefit from the availability of higher-fidelity data sets for a range of systems.

With demand for more calculations of higher precision and increased complexity, comes the need for more complex workflows, 
handled in a systematic, automated manner. To illustrate this point in the context of 
<span style="font-family:american typewriter; font-size:1em;">**exciting**</span>, the choice of the LAPW basis and 
systematic convergence of calculations, even at the ground state level of theory, is more involved than with 
plane wave or Gaussian type orbital (GTO) basis sets. For example, one is free to choose any non-overlapping radii for the muffin-tin spheres of 
each atomic species, and any number of LAPWs and LOs. And for each of these basis functions, one is also free to choose 
the matching order of radial functions, orbital angular momenta, and trial energy parameters associated with them [@gulans:2014]. 
This is before performing the conventional convergence tests, as done in all DFT calculations.

In order to perform systematically-converged, reproducible, benchmark-quality calculations for ground and excited state
phenomena (the latter of which will typically include one or more ground state calculations), a framework is required to assist 
in the selection of calculation parameters, the simplification of input file generation, and the post-processing of 
results. Furthermore, this framework needs to be interoperable with the ecosystem of existing workflow tools. These 
challenges are met by _excitingtools_.

# Summary

_excitingtools_ is a Python package that provides a high-level frontend for the <span style="font-family:american typewriter; font-size:1em;">**exciting**</span>
all-electron DFT package, integrating all aspects of performing a calculation into a single program. _excitingtools_ has
been developed with interoperability in mind, and supports the use of Atomic Simulation Environment (ASE) [@larsen:2017].
Its serialized input classes and output parsers allow it to be used with higher-level workflow managers such as
Jobflow [@Ganose:2022], Atomate [@mathew:2017], and Atomic Simulation Recipes (ASR) [@gjerding:2021]. 


_excitingtools_ is an essential utility for simplifying the use of <span style="font-family:american typewriter; font-size:1em;">**exciting**</span>,
enabling greater user control over calculations. Whilst ASE consists of parsers and calculators, its API is largely 
restricted to ground state energies and forces. _excitingtools_ exposes more functionality allowing users, for instance, 
to analyze their results of different SCF cycles in a calculation, or perform and parse excited state (e.g. GW) calculations. _excitingtools_ enables automation of complex convergence calculations, facilitates high-throughput 
studies, and forms the building blocks of higher-level workflow managers, all of which are prerequisites for moving DFT codes
towards exascale calculations [@gavini:2022]. Furthermore, _excitingtools_ is under active development and follows a continuous 
integration/deployment model, such that new features and updates are delivered several times a year.


## Features

* _excitingtools_ allows the user to quickly create a class object with given key-value pairs in Python, to create input 
  files for <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> in an automated manner.
    - This avoids the need for users to manually configure inputs, which is error-prone, and alleviates the need to 
    frequently write single-purpose scripts.
    - ASE's `Atoms` class is accepted as a structure input.

* _excitingtools_ provides parsers for fifty <span style="font-family:american typewriter; font-size:1em;">**exciting**</span>
    output file formats.
    - Parsing <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> previously required 
    downloading the [NOMAD parsers](https://github.com/nomad-coe/nomad-parser-exciting) which return custom objects, 
    containing copious metadata. This is unnecessary for <span style="font-family:american typewriter; font-size:1em;">**exciting**</span>
    users, and prevents straightforward numerical comparison of parsed results.
    - _excitingtools_ (plus dependencies) is significantly more lightweight than the NOMAD package, approximately two 
    orders of magnitude smaller.
  
* API interoperability and serializable data structures allow easy integration with workflow managers.
    - These features allow the user to create simulation input files, run simulations and analyze data with Python, paving the way to 
    high-throughput calculations with <span style="font-family:american typewriter; font-size:1em;">**exciting**</span>.

# Acknowledgements

The authors would like to thank Ask Hjorth Larsen and David Waroquiers for fruitful discussions. This work received funding from the European Union’s Horizon 2020 research and innovation program under the grant agreement Nº 951786 (NOMAD CoE) and the  German Research Foundation (DFG) through the CRC 1404 (FONDA), project 414984028. Daniel Speckhard acknowledges support by the IMPRS for Elementary Processes in Physical Chemistry.

# References
