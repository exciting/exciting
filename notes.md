## Spin Cases

### spin-unpolarized

* **first variational** (FV) eigenfunctions and eigenvalues

  \[ \psi^{\rm FV}_{n{\bf k}}({\bf r}) = \sum_\mu C^{\rm FV}_{\mu n}({\bf k})\, \phi_{\mu{\bf k}}({\bf r})- [ ] - [x]   \;\text{ and }\; \epsilon^{\rm FV}_{n{\bf k}} \]

* $C^{\rm FV}_{\mu n}({\bf k})$ first variational eigenvectors
  dimensions (per ${\bf k}$-point): basis size $\times$ number of FV states $=$ `ngk+nlotot` $\times$ `nstfv`
* every state $n{\bf k}$ is two-fold degenerate

### spin-polarized

* **second variational** (SV) eigenfunctions and eigenvalues
  \[ \Psi^{\rm SV}_{n{\bf k}}({\bf r}) = \begin{pmatrix} \psi^\uparrow_{n{\bf k}}({\bf r}) \\
  \psi^\downarrow_{n{\bf k}}({\bf r}) \end{pmatrix}
  \;\text{ and }\; \epsilon^{\rm SV}_{n{\bf k}} \]
  $\psi^\sigma_{n{\bf k}}({\bf r})$ with $\sigma = \uparrow / \downarrow$: *spin up* and *spin down* component
  of spinor $\Psi^{\rm SV}_{n{\bf k}}({\bf r})$.
  \[ \psi^\sigma_{n{\bf k}}({\bf r}) = \sum_m C^{\rm SV}_{mn,\sigma}({\bf k})\, \psi^{\rm FV}_{m{\bf k}}({\bf r}) \]
* $C^{\rm SV}_{mn,\sigma}({\bf k})$ second variational eigenvectors
  - dimensions (per ${\bf k}$-point): number of FV states $\times$ number of SV states $\times$ number of spinor components
  $=$ `nstfv` $\times$ `nstsv` $\times$ `nspinor` with `nstsv = nstfv * nspinor` and `nspinor = 2`
  - structure in code: `evecsv(nstsv, nstsv)` for a specific point `ik`
  \[ C^{\rm SV}_{::}({\bf k}) = \begin{pmatrix} C^{\rm SV}_{::,\uparrow}({\bf k}) \\
  C^{\rm SV}_{::,\downarrow}({\bf k}) \end{pmatrix} \]
  where each block is of size `nstfv` $\times$ `nstsv`
* $\epsilon^{\rm SV}_{n{\bf k}}$ second variational eigenvalues
  - dimensions (per ${\bf k}$-point): number of SV states $=$ `nstsv`

#### without spin-orbit coupling / collinear case

* The first `nstfv` states are purely spin up and the second `nstfv` states are purely spin down
  \[ C^{\rm SV}_{::}({\bf k}) = \begin{pmatrix} C^{\rm SV}_{::,\uparrow}({\bf k}) & 0 \\
  0 & C^{\rm SV}_{::,\downarrow}({\bf k}) \end{pmatrix} \]
  where each block is of size `nstfv` $\times$ `nstfv` and
  \[ \epsilon^{\rm SV}_{:{\bf k}} = \begin{pmatrix} \epsilon^\uparrow_{:{\bf k}} \\
  \epsilon^\downarrow_{:{\bf k}} \end{pmatrix} \]
* **Note:** All `nstsv` states are not ordered according to increasing eigenenergy.
  They are separated into two blocks of size `nstfv` for spin up and spin down and 
  ordered according to their eigenenergy within each block.

#### with spin-orbit coupling

* SOC mixes spin components. States are not spin pure but have both spin up and spin down components.
* All `nstsv` states are ordered according to increasing eigenenergy.

## Energy Window Selection
