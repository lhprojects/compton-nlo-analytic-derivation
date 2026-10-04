# Writing conventions

Use these prose, notation, and typesetting rules throughout the book. Physical definitions and normalizations are in `PHYSICS_CONVENTIONS.md`.

- Use sentence case in headings at every level: capitalize only the first word, proper nouns, and acronyms. Keep section headings to one topic. Prefer a prepositional phrase to joining two nouns with "and".
- Short topic headings within a subsection run into the following paragraph. Use `\paragraph{Heading.}` with a final period, and begin the prose on the same line without a forced line or paragraph break after the heading.
- Set group names upright, as in `\mathrm{SU}(2)`. Write the imaginary unit `\ii` (`\mathrm{i}`), but keep an index, variable, or initial-state label `i` italic. Write Hermitian conjugation `\mathrm{h.c.}`.
- Use `D_F` for scalar, `S_F` for Dirac, and `D_{\mu\nu}` (or `D^{\mu\nu}`) for photon Feynman propagators. At first use, state that the photon symbol omits `F`. Use `(S_F)_{ab}` for spinor indices, `S_{F,f}` for species, and `D_{F,0}` for a free scalar with bare mass `m_0`.
- Show the photon gauge parameter as `D_{\mu\nu}^{(\xi)}`. In the Proca comparison, use `D_{\mathrm{cov}}^{\mu\nu}` for the covariant Feynman inverse and `D_T^{\mu\nu}` for the ordinary canonical time-ordered two-point function.
- On first deriving the momentum-space scalar Feynman propagator in Chapter 2, use `\widetilde D_F(p)` to distinguish it from `D_F(x-y)`. After obtaining the expression, state that it is abbreviated `D_F(p)`. Thereafter omit tildes on momentum-space propagators, using the argument to identify the representation. A partial spatial Fourier transform may keep its defined hat.
- Number every displayed formula, including definitions, intermediate steps, boxes, and examples. Use `equation`, `align`, `gather`, etc., never `\[...\]`, starred environments, `\notag`, or `\nonumber`. Number each line of a multiline relation and add `\label` where referenced.
- Avoid semicolons and the words "prescription" and "prescribed" in prose. Name the operation or choice, such as operator ordering, contour choice, or pole displacement. For a source or current, use "given" or "external" when needed, or omit the qualifier. Mathematical notation and code may contain semicolons.
- Avoid "retain" and its forms ("retains", "retained", "retaining") in prose. Use a concrete verb appropriate to the meaning, such as "include", "use", "keep", or "sum".
- Prefer concrete descriptions of where a function vanishes or is nonzero to the mathematical term "support", unless that terminology is needed for the argument.
- When shortening a discussion, remove repeated explanations while preserving the motivation for new terms, assumptions, and definitions. Explain why a construction is useful before introducing it with phrases such as "we take".
- Put standalone diagrams in a numbered `figure` with a caption and text reference. A diagram that must align with an external formula may sit in its equation.
- Lay out Feynman diagrams left to right, with incoming legs left and outgoing legs right. In time-ordered diagrams, time increases rightward and equal-time slices are vertical. Preserve momentum-flow and fermion arrows.
- Put a References section at each chapter's end. Cite numbered references in text, not bibliographic footnotes. Standalone section builds append the shared chapter reference list.
