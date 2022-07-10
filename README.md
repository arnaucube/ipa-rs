# ipa-rs [![Test](https://github.com/arnaucube/ipa-rs/workflows/Test/badge.svg)](https://github.com/arnaucube/ipa-rs/actions?query=workflow%3ATest)

Inner Product Argument (IPA) version from Halo paper (https://eprint.iacr.org/2019/1021.pdf) implementation done to get familiar with [arkworks](https://arkworks.rs) and the modified IPA scheme.


> Warning: do not use this code in production.

### Example

```rust
let mut ipa = IPA::new(8);

let a = vec![
    F::from(1 as u32),
    F::from(2 as u32),
    F::from(3 as u32),
    F::from(4 as u32),
    F::from(5 as u32),
    F::from(6 as u32),
    F::from(7 as u32),
    F::from(8 as u32),
];


let r = F::rand(&mut ipa.rng);

// prover commits
let P = ipa.commit(&a, r).unwrap();


// verifier sets challenges
let U = EdwardsProjective::rand(&mut ipa.rng);
let k = (f64::from(ipa.d as u32).log2()) as usize;
let mut u: Vec<F> = vec![F::zero(); k];
for j in 0..k {
    u[j] = F::rand(&mut ipa.rng);
}
let x = F::from(3 as u32);

// prover opens at the challenges
let b = powers_of(x, ipa.d);
let v = inner_product_field(&a, &b).unwrap();
let proof = ipa.prove(&a, &b, &u, &U).unwrap();

// verifier
let verif = ipa.verify(&x, &v, &P, &proof, &r, &u, &U).unwrap();
assert!(verif);
```
