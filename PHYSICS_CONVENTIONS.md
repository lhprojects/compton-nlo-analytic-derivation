# Physics conventions

Use these definitions and normalizations throughout the book, and define them in the text at first use. Prose and typesetting rules are in `WRITING_STYLE.md`.

## Electric charge

- In four dimensions, $e=|q_e|>0$ is the positive electromagnetic coupling, $q$ (or $q_f=eQ_f$ for species $f$) is signed charge, and $Q_f$ is the dimensionless charge number. Electron and muon fields have $q_e=q_\mu=-e$ and $Q_e=Q_\mu=-1$. Their antiparticles have charge $+e$. Distinguish $Q_f$ from the total charge operator $Q_{\rm em}$.
- Use $D_\mu=\partial_\mu+\mathrm{i}q_f A_\mu$, $\mathcal L_{\rm int}=-q_f\bar\psi_f\gamma^\mu\psi_f A_\mu$, and photon vertex $-\mathrm{i}q_f\gamma^\mu$ ($+\mathrm{i}e\gamma^\mu$ for an electron field). Use the field's $q_f$ for particle and antiparticle external states alike.
- In natural units and the book's electromagnetic normalization, $\alpha=e^2/(4\pi)$ in four dimensions.
- In Chapter 9, use $e$ and $\alpha$ for the physical low-energy charge and fine-structure constant, without a Thomson subscript. Distinguish these from renormalized or running couplings explicitly when those are introduced in later chapters.

## Dimensional regularization

- For dimensional regularization with $d=4-2\epsilon$, renormalized $e$ and $Q_f$ are dimensionless. The interaction coefficient is $\mu^\epsilon eQ_f$, and the full regulated bare coupling is $e_0=\mu^\epsilon Z_e e$, with $q_{0f}=Q_f e_0$. Keep $\mu^\epsilon$ explicit.

## Feynman propagators

- Scalar and Dirac Feynman propagators are vacuum expectations of time-ordered products. Exchanging odd fields in fermionic time ordering contributes a minus sign.
- Propagators include the Feynman-rule factor $\mathrm{i}$. With the mostly-minus metric, the free forms are $D_F(p)=\mathrm{i}/(p^2-m^2+\mathrm{i}0)$, $S_F(p)=\mathrm{i}(\slashed p+m)/(p^2-m^2+\mathrm{i}0)$, and, in Feynman gauge, $D_{\mu\nu}(p)=-\mathrm{i}\eta_{\mu\nu}/(p^2+\mathrm{i}0)$.
- The scalar and spinor subscript $F$ denotes Feynman boundary conditions. The photon propagator has the same boundary conditions despite omitting $F$, while its gauge parameter is a separate choice.
- Relate any auxiliary Dirac matrix without the overall $\mathrm{i}$ to the propagator explicitly. For the free matrix, $S_{0f}(p)=S_{F,f}(p)/\mathrm{i}$, with the Feynman pole displacement understood.
- For a Proca field time-ordered in a chosen inertial frame, the covariant Feynman inverse differs from the ordinary canonical two-point function by a contact term: $D_{\mathrm{cov}}^{\mu\nu}(z)=D_T^{\mu\nu}(z)+\mathrm{i}\delta^\mu{}_0\delta^\nu{}_0\delta^{(4)}(z)/M^2$. The label $\mathrm{cov}$ denotes the covariant kernel.

## Beta functions

- Define $\beta_g=\mu\,\mathrm{d}g/\mathrm{d}\mu$ at fixed bare parameters. For $\ell=\ln(\mu^2/\mu_0^2)$, $\mathrm{d}g/\mathrm{d}\ell=\beta_g/2$. Translate coefficients from other scale-derivative conventions explicitly.
