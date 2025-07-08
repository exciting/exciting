###NOTES

- When you create a macro, or implement a new linalg function, or in general for the whole code
  you need to create two versions one for device and another for host-only version.

- In the compilation, only one set of files depending on the options will be compiled, either ones under 
  device folders (GPU-aware compilations) or the ones under host folders (CPU-aware compilations).

