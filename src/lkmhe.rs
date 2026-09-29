use crate::params::PublicParams;
use qfall_math::integer::{MatZ, Z};
use qfall_math::integer_mod_q::*;
use qfall_math::rational::{MatQ, Q};
use qfall_math::traits::{
    Concatenate, MatrixDimensions, MatrixGetEntry, MatrixGetSubmatrix, MatrixSetEntry,
    MatrixSetSubmatrix, Pow,
};
use qfall_tools::sample::g_trapdoor::{gadget_classical, gadget_parameters::GadgetParameters};

#[derive(Debug, PartialEq, Eq)]
pub enum EvalError {
    MaxLevelReached,
    InvalidLevel,
}

#[derive(Debug, PartialEq, Eq)]
pub enum VectorError {
    WrongShape {
        rows: i64,
        cols: i64,
        expected_cols: usize,
    },
    WrongModulus {
        modulus: Modulus,
        expected_mod: Modulus,
    },
}

#[derive(Clone, Copy)]
pub enum Mode {
    Normal,
    Lossy,
}

pub struct Ciphertext {
    pub c: Vec<MatZq>,
    pub d: MatZq,
    pub t: usize,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SecretKey(MatZq);

impl SecretKey {
    pub fn as_mat(&self) -> &MatZq {
        &self.0
    }

    pub fn add_offset(&self, sigma: &KeyOffset) -> Self {
        Self(&self.0 + &sigma.0)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Message(MatZq);

impl Message {
    pub fn try_from_mat(pp: &PublicParams, mat: MatZq) -> Result<Self, VectorError> {
        validate_vector(pp, &mat)?;
        Ok(Self(mat))
    }

    pub fn sample_uniform(pp: &PublicParams) -> Self {
        Self(MatZq::sample_uniform(1, pp.n, &pp.q))
    }

    pub fn as_mat(&self) -> &MatZq {
        &self.0
    }

    pub fn add_offset(&self, tau: &MessageOffset) -> Self {
        Self(&self.0 + &tau.0)
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct KeyOffset(MatZq);

impl KeyOffset {
    pub fn sample_uniform(pp: &PublicParams) -> Self {
        Self(MatZq::sample_uniform(1, pp.n, &pp.q))
    }

    pub fn as_mat(&self) -> &MatZq {
        &self.0
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MessageOffset(MatZq);

impl MessageOffset {
    pub fn sample_uniform(pp: &PublicParams) -> Self {
        Self(MatZq::sample_uniform(1, pp.n, &pp.q))
    }

    pub fn as_mat(&self) -> &MatZq {
        &self.0
    }
}

/// Returns parameters for toy example. NOT SECURE.
// TODO: implement a function to generate parameters for a given security level
pub fn setup() -> PublicParams {
    let q = Z::from(257);
    let k = 9;
    let m = 16;
    let S = reduced_gadget_basis(m, k, &q);
    let s_gso = MatQ::from(&S).gso();
    PublicParams {
        lambda_1: 8,
        lambda_2: 8,
        t_max: 3,
        q: q,
        n: 2,
        m: m,
        k: k,
        nu: 1.0,
        alpha: 2.0,
        gamma: vec![4.0, 8.0],
        gamma_prime: vec![4.0, 8.0],
        S: S,
        s_gso: s_gso,
    }
}

pub fn keygen(pp: &PublicParams) -> SecretKey {
    SecretKey(MatZq::sample_uniform(1, pp.n, &pp.q))
}

pub fn encrypt(
    pp: &PublicParams,
    s_0: &SecretKey,
    s_1: &SecretKey,
    message: &Message,
    mode: Mode,
) -> Ciphertext {
    let G = reduced_gadget(pp);
    let mut c = Vec::with_capacity(pp.n);
    for i in 0..pp.n {
        let mu_i: Zq = message.as_mat().get_entry(0, i).unwrap();
        let base = build_ct_base(pp, s_0, s_1);
        let ct_i = match mode {
            Mode::Normal => base + &G * mu_i,
            Mode::Lossy => base,
        };
        c.push(ct_i);
    }
    let base = build_ct_base(pp, s_0, s_1);
    let d = match mode {
        Mode::Normal => &base + &G,
        Mode::Lossy => base,
    };
    Ciphertext { c, d, t: 1 }
}

pub fn decrypt(pp: &PublicParams, s_0: &SecretKey, s_1: &SecretKey, ct: &Ciphertext) -> Message {
    let mut mu = MatZq::new(1, pp.n, &pp.q);
    let mut sk = s_0.as_mat().concat_horizontal(s_1.as_mat()).unwrap();
    sk *= -1;
    let one = MatZq::identity(1, 1, &pp.q);
    sk = sk.concat_horizontal(&one).unwrap();
    for i in 0..pp.n {
        let ct_i = &sk * &ct.c[i];
        let v = ct_i.get_submatrix(0, 0, pp.m - pp.k, pp.m - 1).unwrap();
        let mu_i = gadget_decode(pp, &v);
        mu.set_entry(0, i, &mu_i).unwrap();
    }
    Message(mu)
}

pub fn eval(
    pp: &PublicParams,
    sigma_0: &KeyOffset,
    sigma_1: &KeyOffset,
    tau: &MessageOffset,
    ct: &Ciphertext,
) -> Result<Ciphertext, EvalError> {
    if ct.t >= pp.t_max {
        return Err(EvalError::MaxLevelReached);
    } else if ct.t == 0 {
        return Err(EvalError::InvalidLevel);
    }
    let mut c_prime = Vec::with_capacity(pp.n);
    for i in 0..pp.n {
        let c_i = &ct.c[i];
        let mut tau_i: Zq = tau.as_mat().get_entry(0, i).unwrap();
        let mut tau_i_g: MatZq = MatZq::new(1, pp.k, &pp.q);
        for j in 0..pp.k {
            tau_i_g.set_entry(0, j, &tau_i).unwrap();
            tau_i *= 2;
        }
        let mut c_i_prime = c_i + (&ct.d * random_inverse_gadget(pp, &tau_i_g));
        let a_i_prime = c_i_prime.get_submatrix(0, pp.n - 1, 0, pp.m - 1).unwrap();
        let b_i_prime = c_i_prime
            .get_submatrix(pp.n, 2 * pp.n - 1, 0, pp.m - 1)
            .unwrap();
        let sigma_0_a_i_prime = sigma_0.as_mat() * a_i_prime;
        let sigma_1_b_i_prime = sigma_1.as_mat() * b_i_prime;
        let e_i = MatZ::sample_discrete_gauss(1, pp.m, 0, pp.gamma[ct.t - 1]).unwrap();
        let mut bottom = sigma_0_a_i_prime + sigma_1_b_i_prime + e_i;
        bottom += c_i_prime.get_row(2 * pp.n).unwrap();
        c_i_prime.set_row(2 * pp.n, &bottom, 0).unwrap();
        c_prime.push(c_i_prime);
    }
    let g = gadget_classical::gen_gadget_vec(pp.k, &Z::from(2)).transpose();
    let g_mod_q = MatZq::from((&g, &pp.q));
    let mut d_prime = &ct.d * random_inverse_gadget(pp, &g_mod_q);
    let a_prime = d_prime.get_submatrix(0, pp.n - 1, 0, pp.m - 1).unwrap();
    let b_prime = d_prime
        .get_submatrix(pp.n, 2 * pp.n - 1, 0, pp.m - 1)
        .unwrap();
    let sigma_0_a_prime = sigma_0.as_mat() * a_prime;
    let sigma_1_b_prime = sigma_1.as_mat() * b_prime;
    let e = MatZ::sample_discrete_gauss(1, pp.m, 0, pp.gamma_prime[ct.t - 1]).unwrap();
    let mut bottom = sigma_0_a_prime + sigma_1_b_prime + e;
    bottom += d_prime.get_row(2 * pp.n).unwrap();
    d_prime.set_row(2 * pp.n, &bottom, 0).unwrap();

    Ok(Ciphertext {
        c: c_prime,
        d: d_prime,
        t: ct.t + 1,
    })
}

fn validate_vector(pp: &PublicParams, value: &MatZq) -> Result<(), VectorError> {
    let rows = value.get_num_rows();
    let cols = value.get_num_columns();
    if rows != 1 || cols != i64::try_from(pp.n).unwrap() {
        return Err(VectorError::WrongShape {
            rows,
            cols,
            expected_cols: pp.n,
        });
    } else if value.get_mod() != Modulus::from(&pp.q) {
        return Err(VectorError::WrongModulus {
            modulus: value.get_mod(),
            expected_mod: Modulus::from(&pp.q),
        });
    }
    Ok(())
}

fn reduced_gadget(pp: &PublicParams) -> MatZq {
    let mut G = MatZq::new(2 * pp.n + 1, pp.m, &pp.q);
    let mut value = Zq::from((1, &pp.q));
    for i in (pp.m - pp.k)..pp.m {
        G.set_entry(2 * pp.n, i, &value).unwrap();
        value *= 2;
    }
    G
}

fn build_ct_base(pp: &PublicParams, s_0: &SecretKey, s_1: &SecretKey) -> MatZq {
    // sample error
    let e_i = MatZ::sample_discrete_gauss(1, pp.m, 0, pp.nu).unwrap();
    let a_i = MatZq::sample_uniform(pp.n, pp.m, &pp.q);
    let b_i = MatZq::sample_uniform(pp.n, pp.m, &pp.q);
    let ab_i = a_i.concat_vertical(&b_i).unwrap();
    let bottom = s_0.as_mat() * &a_i + s_1.as_mat() * &b_i + &e_i;
    ab_i.concat_vertical(&bottom).unwrap()
}

fn reduced_gadget_basis(m: usize, k: usize, q: &Z) -> MatZ {
    let params = GadgetParameters::init_default(1, q);
    let S_k = gadget_classical::short_basis_gadget(&params);
    let unit = MatZ::identity(m - k, m - k);
    let zero = MatZ::new(m - k, k);
    let bottom = MatZ::new(k, m - k).concat_horizontal(&S_k).unwrap();
    let S = unit.concat_horizontal(&zero).unwrap();
    S.concat_vertical(&bottom).unwrap()
}

fn bit_decomp(pp: &PublicParams, a: &Zq) -> MatZ {
    // q must be smaller than 64 bits
    // TODO: rewrite to avoid conversion into u64
    let a: u64 = a
        .get_representative_least_nonnegative_residue()
        .try_into()
        .unwrap();
    let mut bits = MatZ::new(1, pp.k);
    for i in 0..pp.k {
        let bit = (a >> i) & 1;
        bits.set_entry(0, i, &Z::from(bit)).unwrap();
    }
    bits
}

fn random_inverse_gadget(pp: &PublicParams, v: &MatZq) -> MatZ {
    let mut X = MatZ::new(pp.m, pp.m);
    for i in 0..pp.m - pp.k {
        let c_i = MatZ::sample_d_precomputed_gso(&pp.S, &pp.s_gso, &MatQ::new(pp.m, 1), pp.alpha)
            .unwrap();
        X.set_column(i, &c_i, 0).unwrap();
    }
    for i in pp.m - pp.k..pp.m {
        let j = i - (pp.m - pp.k);
        let v_i: Zq = v.get_entry(0, j).unwrap();
        let bits = bit_decomp(pp, &v_i);
        let mut c_0 = MatZ::new(pp.m, 1);
        for l in 0..pp.k {
            c_0.set_entry(pp.m - pp.k + l, 0, &bits.get_entry(0, l).unwrap())
                .unwrap();
        }
        let center: MatQ = c_0.clone().try_into().unwrap();
        let center = Z::from(-1) * center;
        let y = MatZ::sample_d_precomputed_gso(&pp.S, &pp.s_gso, &center, pp.alpha).unwrap();
        let c_i = y + c_0;
        X.set_column(i, &c_i, 0).unwrap();
    }
    X
}

fn gadget_decode(pp: &PublicParams, v: &MatZq) -> Zq {
    let mu: Zq;
    let mut T = Z::ZERO;
    for i in 0..&pp.k - 1 {
        let v_i: Zq = v.get_entry(0, i).unwrap();
        let v_next: Zq = v.get_entry(0, i + 1).unwrap();
        let delta: Zq = 2 * v_i - v_next;
        // let d = delta.get_representative_least_absolute_residue();
        // there is a bug in get_representative_least_absolute_residue() so compute centered representative manually
        let mut d = delta.get_representative_least_nonnegative_residue();
        if &d * 2 > pp.q {
            d -= &pp.q;
        }
        T = 2 * T + d;
    }
    let denominator = Z::from(2).pow(u64::try_from(pp.k - 1).unwrap()).unwrap();

    let e_0 = (Q::from(T) / Q::from(denominator)).round();
    let v_0: Zq = v.get_entry(0, 0).unwrap();
    mu = v_0 - Zq::from((e_0, &pp.q));
    mu
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gadget_basis() {
        let pp = setup();
        let S = reduced_gadget_basis(pp.m, pp.k, &pp.q);
        let G = reduced_gadget(&pp);
        let zero = MatZq::new(2 * pp.n + 1, pp.m, &pp.q);
        assert_eq!(G * &S, zero);
        let det = S.det().unwrap();
        assert_eq!(det.abs(), pp.q);
    }

    #[test]
    fn bit_decomp_test() {
        let pp = setup();
        let v = Zq::sample_uniform(&pp.q);
        let bits = bit_decomp(&pp, &v);
        let mut c_0 = MatZ::new(pp.m, 1);
        for l in 0..pp.k {
            c_0.set_entry(pp.m - pp.k + l, 0, &bits.get_entry(0, l).unwrap())
                .unwrap();
        }
        let G = reduced_gadget(&pp);
        let mut G_c_0 = MatZq::new(2 * pp.n + 1, 1, &pp.q);
        G_c_0.set_entry(2 * pp.n, 0, &v).unwrap();
        assert_eq!(G * c_0, G_c_0);
    }

    #[test]
    fn random_inverse_gadget_test() {
        let pp = setup();
        let v = MatZq::sample_uniform(1, pp.k, &pp.q);
        let X = random_inverse_gadget(&pp, &v);
        let G = reduced_gadget(&pp);
        let mut G_X = MatZq::new(2 * pp.n + 1, pp.m, &pp.q);
        for i in 0..pp.k {
            let v_i: Zq = v.get_entry(0, i).unwrap();
            G_X.set_entry(2 * pp.n, pp.m - pp.k + i, &v_i).unwrap();
        }
        assert_eq!(G * X, G_X);
    }

    #[test]
    fn correctness() {
        let pp = setup();
        let s_0 = keygen(&pp);
        let s_1 = keygen(&pp);
        let message = Message::sample_uniform(&pp);
        let ct = encrypt(&pp, &s_0, &s_1, &message, Mode::Normal);
        let decrypted = decrypt(&pp, &s_0, &s_1, &ct);
        assert_eq!(message, decrypted);
    }

    #[test]
    fn kmhe_correctness_normal() {
        let pp = setup();
        let s_0 = keygen(&pp);
        let s_1 = keygen(&pp);
        let sigma_0 = KeyOffset::sample_uniform(&pp);
        let sigma_1 = KeyOffset::sample_uniform(&pp);
        let tau = MessageOffset::sample_uniform(&pp);
        let m = Message::sample_uniform(&pp);
        let ct_1 = encrypt(&pp, &s_0, &s_1, &m, Mode::Normal);
        let ct_2 = eval(&pp, &sigma_0, &sigma_1, &tau, &ct_1).unwrap();
        let s_0_prime = s_0.add_offset(&sigma_0);
        let s_1_prime = s_1.add_offset(&sigma_1);
        let m_prime = m.add_offset(&tau);
        let m_star = decrypt(&pp, &s_0_prime, &s_1_prime, &ct_2);
        assert_eq!(m_star, m_prime);
    }

    #[test]
    fn kmhe_correctness_lossy() {
        let pp = setup();
        let s_0 = keygen(&pp);
        let s_1 = keygen(&pp);
        let sigma_0 = KeyOffset::sample_uniform(&pp);
        let sigma_1 = KeyOffset::sample_uniform(&pp);
        let tau = MessageOffset::sample_uniform(&pp);
        let m = Message::sample_uniform(&pp);
        let ct_1 = encrypt(&pp, &s_0, &s_1, &m, Mode::Lossy);
        let ct_2 = eval(&pp, &sigma_0, &sigma_1, &tau, &ct_1).unwrap();
        let s_0_prime = s_0.add_offset(&sigma_0);
        let s_1_prime = s_1.add_offset(&sigma_1);
        let m_prime = Message::try_from_mat(&pp, MatZq::new(1, pp.n, &pp.q)).unwrap();
        let m_star = decrypt(&pp, &s_0_prime, &s_1_prime, &ct_2);
        assert_eq!(m_star, m_prime);
    }
}
