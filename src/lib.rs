extern crate ark_ed_on_bn254;
use ark_ec::ProjectiveCurve;
use ark_ed_on_bn254::{EdwardsProjective, Fr};
use ark_ff::{fields::PrimeField, Field};
use ark_std::{One, UniformRand, Zero};

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

    pub fn commit(&self, a: &[Fr], r: Fr) -> Result<EdwardsProjective, String> {
        Ok(inner_product_point(a, &self.Gs)? + self.H.mul(r.into_repr()))
    }

    pub fn ipa(
        &mut self,
        a: &[Fr],
        b: &[Fr],
        u: &[Fr],
        U: &EdwardsProjective,
    ) -> Result<Proof, String> {
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

            L[j] = inner_product_point(&a_lo, &G_hi)?
                + self.H.mul(l[j].into_repr())
                + U.mul(inner_product_field(&a_lo, &b_hi)?.into_repr());
            R[j] = inner_product_point(&a_hi, &G_lo)?
                + self.H.mul(r[j].into_repr())
                + U.mul(inner_product_field(&a_hi, &b_lo)?.into_repr());

            let uj = u[j];
            let uj_inv = u[j].inverse().unwrap();

            a = vec_add(
                &vec_scalar_mul_field(&a_lo, &uj),
                &vec_scalar_mul_field(&a_hi, &uj_inv),
            )?;
            b = vec_add(
                &vec_scalar_mul_field(&b_lo, &uj_inv),
                &vec_scalar_mul_field(&b_hi, &uj),
            )?;
            G = vec_add_point(
                &vec_scalar_mul_point(&G_lo, &uj_inv),
                &vec_scalar_mul_point(&G_hi, &uj),
            )?;
        }

        if a.len() != 1 {
            return Err(format!("a.len() should be 1, a.len()={}", a.len()));
        }
        if b.len() != 1 {
            return Err(format!("b.len() should be 1, b.len()={}", b.len()));
        }
        if G.len() != 1 {
            return Err(format!("G.len() should be 1, G.len()={}", G.len()));
        }

        Ok(Proof {
            a: a[0],
            l,
            r,
            L,
            R,
        })
    }
    pub fn verify(
        &self,
        x: &Fr,
        P: &EdwardsProjective,
        p: &Proof,
        r: &Fr,
        u: &[Fr],
        U: &EdwardsProjective,
    ) -> Result<bool, String> {
        let mut q_0 = *P;
        let mut r = *r;

        // compute b & G from s
        let s = build_s(u, self.d as usize);
        let bs = powers_of(*x, self.d);
        let b = inner_product_field(&s, &bs)?;
        let G = inner_product_point(&s, &self.Gs)?;

        #[allow(clippy::needless_range_loop)]
        for j in 0..u.len() {
            let uj2 = u[j].square();
            let uj_inv2 = u[j].inverse().unwrap().square();

            q_0 = q_0 + p.L[j].mul(uj2.into_repr()) + p.R[j].mul(uj_inv2.into_repr());
            r = r + p.l[j] * uj2 + p.r[j] * uj_inv2;
        }

        let q_1 = G.mul(p.a.into_repr()) + self.H.mul(r.into_repr()) + U.mul((p.a * b).into_repr());

        Ok(q_0 == q_1)
    }
}

// s = (
//   u₁⁻¹ u₂⁻¹ … uₖ⁻¹,
//   u₁   u₂⁻¹ … uₖ⁻¹,
//   u₁⁻¹ u₂   … uₖ⁻¹,
//   u₁   u₂   … uₖ⁻¹,
//   ⋮    ⋮      ⋮
//   u₁   u₂   … uₖ
// )
fn build_s(u: &[Fr], d: usize) -> Vec<Fr> {
    let k = (f64::from(d as u32).log2()) as usize;
    let mut s: Vec<Fr> = vec![Fr::one(); d];
    let mut t = d;
    for j in (0..k).rev() {
        t /= 2;
        let mut c = 0;
        #[allow(clippy::needless_range_loop)]
        for i in 0..d {
            if c < t {
                s[i] *= u[j].inverse().unwrap();
            } else {
                s[i] *= u[j];
            }
            c += 1;
            if c >= t * 2 {
                c = 0;
            }
        }
    }
    s
}

fn inner_product_field(a: &[Fr], b: &[Fr]) -> Result<Fr, String> {
    if a.len() != b.len() {
        return Err(format!(
            "a.len()={} must be equal to b.len()={}",
            a.len(),
            b.len()
        ));
    }
    let mut c: Fr = Fr::zero();
    for i in 0..a.len() {
        c += a[i] * b[i];
    }
    Ok(c)
}

fn inner_product_point(a: &[Fr], b: &[EdwardsProjective]) -> Result<EdwardsProjective, String> {
    if a.len() != b.len() {
        return Err(format!(
            "a.len()={} must be equal to b.len()={}",
            a.len(),
            b.len()
        ));
    }
    let mut c: EdwardsProjective = EdwardsProjective::zero();
    for i in 0..a.len() {
        c += b[i].mul(a[i].into_repr());
    }
    Ok(c)
}

fn vec_add(a: &[Fr], b: &[Fr]) -> Result<Vec<Fr>, String> {
    if a.len() != b.len() {
        return Err(format!(
            "a.len()={} must be equal to b.len()={}",
            a.len(),
            b.len()
        ));
    }
    let mut c: Vec<Fr> = vec![Fr::zero(); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] + b[i];
    }
    Ok(c)
}
fn vec_add_point(
    a: &[EdwardsProjective],
    b: &[EdwardsProjective],
) -> Result<Vec<EdwardsProjective>, String> {
    if a.len() != b.len() {
        return Err(format!(
            "a.len()={} must be equal to b.len()={}",
            a.len(),
            b.len()
        ));
    }
    let mut c: Vec<EdwardsProjective> = vec![EdwardsProjective::zero(); a.len()];
    for i in 0..a.len() {
        c[i] = a[i] + b[i];
    }
    Ok(c)
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

fn powers_of(x: Fr, d: u32) -> Vec<Fr> {
    let mut c: Vec<Fr> = vec![Fr::zero(); d as usize];
    c[0] = x;
    for i in 1..d as usize {
        c[i] = c[i - 1] * x;
    }
    c
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;

    #[test]
    fn test_utils() {
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
        let c = inner_product_field(&a, &b).unwrap();
        assert_eq!(c, Fr::from(30 as u32));
    }

    #[test]
    fn test_homomorphic_property() {
        let d = 8;
        let ipa = IPA::new(d);

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
        let b = a.clone();

        let mut rng = ark_std::rand::thread_rng();
        let r = Fr::rand(&mut rng);
        let s = Fr::rand(&mut rng);

        let vc_a = ipa.commit(&a, r).unwrap();
        let vc_b = ipa.commit(&b, s).unwrap();

        let expected_vc_c = ipa.commit(&vec_add(&a, &b).unwrap(), r + s).unwrap();
        let vc_c = vc_a + vc_b;
        assert_eq!(vc_c, expected_vc_c);
    }

    #[test]
    fn test_inner_product_argument_proof() {
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

        let mut P = ipa.commit(&a, r).unwrap();
        let v = inner_product_field(&a, &b).unwrap();

        let U = EdwardsProjective::rand(&mut ipa.rng);

        let k = (f64::from(ipa.d as u32).log2()) as usize;
        let mut u: Vec<Fr> = vec![Fr::zero(); k];
        for j in 0..k {
            u[j] = Fr::rand(&mut ipa.rng);
        }

        P = P + U.mul(v.into_repr());

        let proof = ipa.ipa(&a, &b, &u, &U).unwrap();
        let verif = ipa.verify(&x, &P, &proof, &r, &u, &U).unwrap();
        assert!(verif);
    }
}
