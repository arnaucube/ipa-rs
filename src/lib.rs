extern crate ark_ed_on_bn254;
use ark_ec::ProjectiveCurve;
use ark_ed_on_bn254::{EdwardsProjective, Fr};
use ark_ff::{fields::PrimeField, Field}; // BigInteger
use ark_std::{UniformRand, Zero};

#[allow(non_snake_case)]
pub struct IPA {
    d: u32,
    H: EdwardsProjective,
    Gs: Vec<EdwardsProjective>,
    rng: rand::rngs::ThreadRng,
}

#[allow(non_snake_case)]
pub struct Proof {
    a: Fr,
    b: Fr,                // TODO not needed
    G: EdwardsProjective, // TODO not needed
    l: Vec<Fr>,
    r: Vec<Fr>,
    L: Vec<EdwardsProjective>,
    R: Vec<EdwardsProjective>,
}

#[allow(non_snake_case)]
#[allow(clippy::many_single_char_names)]
impl IPA {
    pub fn new(d: u32) -> IPA {
        let mut rng = ark_std::rand::thread_rng();

        let mut gs: Vec<EdwardsProjective> = Vec::new();
        for _ in 0..d {
            gs.push(EdwardsProjective::rand(&mut rng));
        }

        IPA {
            d,
            H: EdwardsProjective::rand(&mut rng),
            Gs: gs,
            rng,
        }
    }

    pub fn commit(&self, a: &[Fr], r: Fr) -> EdwardsProjective {
        inner_product_point(a, &self.Gs) + self.H.mul(r.into_repr())
    }

    pub fn ipa(&mut self, a: &[Fr], b: &[Fr], u: &[Fr], U: &EdwardsProjective) -> Proof {
        let mut a = a.to_owned();
        let mut b = b.to_owned();
        let mut G = self.Gs.clone();

        let k = (f64::from(self.d as u32).log2()) as usize;
        let mut l: Vec<Fr> = vec![Fr::zero(); k];
        let mut r: Vec<Fr> = vec![Fr::zero(); k];
        let mut L: Vec<EdwardsProjective> = vec![EdwardsProjective::zero(); k];
        let mut R: Vec<EdwardsProjective> = vec![EdwardsProjective::zero(); k];

        for j in (0..k).rev() {
            let m = a.len() / 2;
            let a_lo = a[..m].to_vec();
            let a_hi = a[m..].to_vec();
            let b_lo = b[..m].to_vec();
            let b_hi = b[m..].to_vec();
            let G_lo = G[..m].to_vec();
            let G_hi = G[m..].to_vec();

            l[j] = Fr::rand(&mut self.rng);
            r[j] = Fr::rand(&mut self.rng);

            L[j] = inner_product_point(&a_lo, &G_hi)
                + self.H.mul(l[j].into_repr())
                + U.mul(inner_product_field(&a_lo, &b_hi).into_repr());
            R[j] = inner_product_point(&a_hi, &G_lo)
                + self.H.mul(r[j].into_repr())
                + U.mul(inner_product_field(&a_hi, &b_lo).into_repr());

            let uj = u[j];
            let uj_inv = u[j].inverse().unwrap();

            a = vec_add(
                &vec_scalar_mul_field(&a_lo, &uj),
                &vec_scalar_mul_field(&a_hi, &uj_inv),
            );
            b = vec_add(
                &vec_scalar_mul_field(&b_lo, &uj_inv),
                &vec_scalar_mul_field(&b_hi, &uj),
            );
            G = vec_add_point(
                &vec_scalar_mul_point(&G_lo, &uj_inv),
                &vec_scalar_mul_point(&G_hi, &uj),
            );
        }

        // TODO assert len a,b,G == 1

        Proof {
            a: a[0],
            b: b[0],
            G: G[0],
            l,
            r,
            L,
            R,
        }
    }
    pub fn verify(
        &self,
        P: &EdwardsProjective,
        p: &Proof,
        r: &Fr,
        u: &[Fr],
        U: &EdwardsProjective,
    ) -> bool {
        let mut q_0 = *P;
        let mut r = *r;

        // TODO compute b & G without getting them in the proof package

        #[allow(clippy::needless_range_loop)]
        for j in 0..u.len() {
            let uj2 = u[j].square();
            let uj_inv2 = u[j].inverse().unwrap().square();

            q_0 = q_0 + p.L[j].mul(uj2.into_repr()) + p.R[j].mul(uj_inv2.into_repr());
            r = r + p.l[j] * uj2 + p.r[j] * uj_inv2;
        }

        let q_1 =
            p.G.mul(p.a.into_repr()) + self.H.mul(r.into_repr()) + U.mul((p.a * p.b).into_repr());

        q_0 == q_1
    }
}

fn inner_product_field(a: &[Fr], b: &[Fr]) -> Fr {
    // TODO require lens equal
    let mut c: Fr = Fr::zero();
    for i in 0..a.len() {
        c += a[i] * b[i];
    }
    c
}

fn inner_product_point(a: &[Fr], b: &[EdwardsProjective]) -> EdwardsProjective {
    // TODO require lens equal
    let mut c: EdwardsProjective = EdwardsProjective::zero();
    for i in 0..a.len() {
        c += b[i].mul(a[i].into_repr());
    }
    c
}

fn vec_add(a: &[Fr], b: &[Fr]) -> Vec<Fr> {
    // TODO require len equal
    let mut c: Vec<Fr> = vec![Fr::zero(); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] + b[i];
    }
    c
}
fn vec_add_point(a: &[EdwardsProjective], b: &[EdwardsProjective]) -> Vec<EdwardsProjective> {
    // TODO require len equal
    let mut c: Vec<EdwardsProjective> = vec![EdwardsProjective::zero(); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] + b[i];
    }
    c
}

fn vec_scalar_mul_field(a: &[Fr], b: &Fr) -> Vec<Fr> {
    let mut c: Vec<Fr> = vec![Fr::zero(); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] * b;
    }
    c
}
fn vec_scalar_mul_point(a: &[EdwardsProjective], b: &Fr) -> Vec<EdwardsProjective> {
    let mut c: Vec<EdwardsProjective> = vec![EdwardsProjective::zero(); a.len()];
    for i in 0..a.len() {
        c[i] = a[i].mul(b.into_repr());
    }
    c
}

#[allow(dead_code)]
fn powers_of(x: Fr, d: u32) -> Vec<Fr> {
    let mut c: Vec<Fr> = vec![Fr::zero(); d as usize];
    c[0] = x;
    for i in 1..d as usize {
        // TODO redo better
        c[i] = c[i - 1] * x;
    }
    c
}

// fn inner_product<T>(a: Vec<T>, b: Vec<T>) -> T {
//     // require lens equal
//     let mut c: T = Zero();
//     for i in 0..a.len() {
//         c = c + a[i] * b[i];
//     }
//     c
// }

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;

    #[test]
    fn test_utils() {
        // let a = Fr::from(1 as u32);
        // let b = Fr::one();
        // println!("A: {:?}", Fr::from(1 as u32));
        // println!("A: {:?}", a);
        // println!("B: {:?}", b);

        let a = vec![
            Fr::from(1 as u32),
            Fr::from(2 as u32),
            Fr::from(3 as u32),
            Fr::from(4 as u32),
        ];
        let b = vec![
            Fr::from(1 as u32),
            Fr::from(2 as u32),
            Fr::from(3 as u32),
            Fr::from(4 as u32),
        ];
        let c = inner_product_field(&a, &b);
        println!("c: {:?}", c);

        // let result = 2 + 2;
        // assert_eq!(result, 4);
    }

    #[test]
    fn test_inner_product() {
        let d = 8;
        let mut ipa = IPA::new(d);

        let a = vec![
            Fr::from(1 as u32),
            Fr::from(2 as u32),
            Fr::from(3 as u32),
            Fr::from(4 as u32),
            Fr::from(5 as u32),
            Fr::from(6 as u32),
            Fr::from(7 as u32),
            Fr::from(8 as u32),
        ];

        let x = Fr::from(3 as u32);
        let b = powers_of(x, ipa.d);

        let r = Fr::rand(&mut ipa.rng);

        let mut P = ipa.commit(&a, r);
        let v = inner_product_field(&a, &b);

        let U = EdwardsProjective::rand(&mut ipa.rng);

        let k = (f64::from(ipa.d as u32).log2()) as usize;
        let mut u: Vec<Fr> = vec![Fr::zero(); k];
        for j in 0..k {
            u[j] = Fr::rand(&mut ipa.rng);
        }

        P = P + U.mul(v.into_repr());

        let proof = ipa.ipa(&a, &b, &u, &U);
        let verif = ipa.verify(&P, &proof, &r, &u, &U);
        assert!(verif);
    }
}
