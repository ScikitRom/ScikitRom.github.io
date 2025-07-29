### **scikit-rom**

*scikit-rom* is a lightweight, Python-based platform designed for projection-based model reduction of finite element models with moderate to large problem sizes (up to approximately 600,000 degrees of freedom). Built on top of `scikit-fem`, an easy-to-use finite element assembly library, this repository demonstrates reduced-order modeling (ROM) techniques for a range of linear and nonlinear problems, including time-dependent ones. Current examples primarily focus on thermal and mechanical systems, with extensions to fluid dynamics planned for future releases.

Beyond its application to real-world finite element problems, the library also serves as a valuable testing ground for new model reduction algorithms and as a tutorial resource for those new to the field.

---

### **What is Reduced-Order Modeling?**

Reduced-order modeling (ROM) is a computational strategy for approximating high-fidelity numerical models with significantly reduced computational cost. Instead of solving large-scale systems—often with millions of unknowns—ROM techniques construct low-dimensional surrogates that retain the dominant physics of the system, enabling efficient and repeated evaluations.

The underlying principle is that many physical systems exhibit behavior that can be captured by a small number of dominant modes, despite their discretization on fine spatial meshes. Techniques such as Proper Orthogonal Decomposition (POD) extract these modes from full-order simulations, enabling the construction of reduced models that are both accurate and computationally efficient.

---

### **Offline–Online Decomposition**

*scikit-rom* implements the classical offline–online decomposition paradigm:

* **Offline Phase** (computationally intensive, performed once):

  * Sampling of the parameter space
  * Generation of full-order solution snapshots
  * Computation of the reduced basis via singular value decomposition (SVD)
  * Assembly of parameter-independent reduced operators

* **Online Phase** (efficient, performed repeatedly):

  * Evaluation of the system response for new parameter values using preassembled reduced operators
  * Speedups on the order of 10–500× relative to full-order models are typical