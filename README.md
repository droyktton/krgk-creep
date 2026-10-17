# creep-T0plus-bubble

T → 0⁺ creep dynamics (Kolton, Rosso, Giamarchi, Krauth, PRL 97, 057001 (2006); PRB 79, 184207 (2009); Ferrero et al., PRL 118, 147208 (2017)) for a straight wall and for an expanding circular domain (bubble), with the same energy and the same isotropic disorder, to test whether the bubble geometry changes the roughness exponent in deep creep (H ≪ H_d), where Langevin dynamics cannot go.

## Motivation
Expanding bubbles in creep give ζ ≈ 0.6 (Durin et al., PRB 110, L020405 (2024); our CoNi bubble), straight creep walls give ζ ≈ 1.0–1.25 (driven-magnetic-DW-persistence, UNIVERSALITY.md). In φ⁴ (T = 0 near h_c, and Langevin at T > 0 down to h/h_c ≈ 0.6) bubble and stripe of the same age give the same ζ. This code reaches H/H_d → 0.

## Contents
- `original/vmc/`: the KRGK code as found on the cluster (`~/Codigos/cuerdas_elasticas/vmc`, bitbucket `abk_balde/cuerdas_elasticas`, last change 2019). C++, Blitz++, serial. Not modified.
- `original/common/`: shared helpers of that repository (the large safe-prime tables were left out).

## Plan
1. Build and run the original (straight line, harmonic elasticity, lattice disorder) and reproduce L_opt(F) and the creep events of KRGK/Ferrero as a reference.
2. Write the new code with a geometry-agnostic core (the energy-ordered search of `Besc_new2` and the forward relaxation), and two geometries:
   - straight: u_i on x_i = i·a, periodic in x;
   - bubble: R_i on θ_i = 2πi/N, periodic in θ.
3. Energy, same form in both geometries (arc element a_i = a or R_i·Δθ):
   - elastic: c·Σ √(a_i² + (Δ_i)²), or its harmonic expansion c·Σ [a_i + Δ_i²/(2a_i)]; the bubble gets the Laplace pressure from the first term;
   - drive: −H·(swept area): −H·a·Σ u_i, or −H·Δθ·Σ R_i²/2;
   - pinning: Σ a_i·V(r_i), with V an isotropic random field defined in the plane (correlation length ξ), sampled at r_i = (x_i, u_i) or (R_i cos θ_i, R_i sin θ_i). The disorder does not live on the computational grid, so the polar grid adds no anisotropy.
4. Discretization: transverse step δ ≪ ξ (u or R in units of δ, |Δ_i| ≤ Mcut as in the original); arc spacing a_i ≲ ξ, with N fixed by the final radius, or doubled as R grows.
5. Measurements: w²(ℓ) and ζ at equal age, P_ov, events (L_opt, avalanche sizes, cluster shapes S ~ ℓ^(1+ζ)), with the tools of driven-magnetic-DW-persistence.

No velocity enters the T → 0⁺ dynamics, so there is no geometric KPZ term: the model is qEW-like in both geometries, and a bubble/straight difference would be purely geometric (dilution, curvature).
