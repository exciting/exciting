In general, the crystal under consideration exhibits certain symmetry
operations \(\mathcal{S}\) under which the density and potential are
invariant. All symmetry operations form a group called \(\mathcal{G}\).
Each symmetry operation \(\mathcal{S}=\lbrace {\bf \rm S}, {\bf \tau}_S \rbrace\)
is formed by an (im)proper rotation matrix \({\bf \rm S}\) and a translation
vector \({\bf \tau}_S\) such that each atom \(\alpha\) is mapped into an 
chemically equivalent atom \(\alpha'\), i.e,
\[ \mathcal{S}\, {\bf \tau}_\alpha = {\bf \rm S} ({\bf \tau}_\alpha + {\bf \tau}_S)
   = {\bf \tau}_{\alpha'} + {\bf R}_S \;.\]
@note Mostly, in the literature, the symmetry operations are defined by the rotation
followed by the translation. In \(\texttt{exciting}\), they are defined by the translation
followed by the rotation. @endnote

The goal is to exploit the symmetries of the crystal to reduce the
computational costs. To this end, we want to reduce the number of 
\({\bf k}\) points in the construction of the density response to the 
*irreducible Brillouin zone* (IBZ), i.e., the set of \({\bf k}\) points
from which all other points in the Brillouin zone (BZ) can be reached
by rotations with the crystal symmetries. The missing information can then
be restored by symmetrizing the density response with all crystal symmetries.

Let's consider the perturbation of a phonon-like displacement of atom \(\alpha\) 
in the Cartesian direction \(i\) with wavevector \({\bf q}\), 
\(\delta^{\bf q}_{\alpha i}\). The problem is that the application of each
symmetry operation mixes phonon wavectors, atoms and Cartesian directions
\[ \mathcal{S}\, \delta^{\bf q}_{\alpha i} 
   \rightarrow \sum_{j=1}^3 S_{ij}\, \delta^{{\bf q}'}_{\alpha' j} \;.\]
As a consequence, in order to symmetrize the density response, one has to calculate
the response to \(N_{\bf q} \times 3 N_{\rm at}\) simultaneously which is 
impractial in terms of memory consumption. 

The first step to overcome this problem
is to reduce the group of symmetry operations to the ones that leave \({\bf q}\) invariant,
the so-called *small group of* \({\bf q}\)
\[ \mathcal{G}_{\bf q} = \lbrace \mathcal{S} \in \mathcal{G} : {\bf q}' = {\bf \rm S}{\bf q}
   = {\bf q} + {\bf G} \rbrace \;. \]
This avoids the mixing of \({\bf q}\) vectors and reduces the number of simultaneous 
perturbations to \(3N_{\rm at}\). The cost for this gain is a generally increased number
of \({\bf k}\) vectors to consider in the construction of the density response.
Due to the reduced number of symmetry operations, \({\bf k}\) can only be reduced
to a \({\bf q}\)-dependent increased irreducible Brillouin zone \(\text{IBZ}({\bf q})\).

The second step is to use special symmetry adapted displacement patterns \(\delta^{\rm q}_{I \mu}\)
expressed as linear combinations of the canonical displacement patterns
\[ \delta^{\bf q}_{I \mu} = \sum_{\alpha,i} p^{I \mu}_{\alpha i}({\bf q})\, \delta^{\bf q}_{\alpha i} \;. \]
The coefficients \(p^{I \mu}_{\alpha i}({\bf q})\) describe the so-called 
*irreducible representations* (irreps) \(I\). The index \(\mu = 1, \dots, d_I\) labels
the *members of the irrep* \(I\) and \(d_I\) is the *dimension of the irrep* \(I\).
The dimension is typically a relatively small number \(\leq 3\) and at most it is 6.
Now, the symmetry operations only mix the members of an irrep
\[ \mathcal{S}\, \delta^{\bf q}_{I \mu} 
   = \sum_{\nu=1}^{d_I} \texttt{S}^I_{\mu \nu}({\bf q})\, \delta^{\bf q}_{I \nu} \;, \]
where \({\bf \texttt{S}}^I({\bf q})\) is the matrix representation of the symmetry operation
\(\mathcal{S}\) in the basis of the irrep \(I\). 

Using the irrep displacement patterns,
each pair \(({\bf q},I)\) forms an independent part of the full phonon calculation
and within each part, only \(d_I\) perturbations have to be treated simultaneously.

An extensive review of symmetries in phonons can be found in
[A. Maradudin and S. Vosko. Symmetry Properties of the Normal Vibrations of a Crystal, 
*Rev. Mod. Phys.* **40**, 1 (1968)](https://doi.org/10.1103/RevModPhys.40.1)
