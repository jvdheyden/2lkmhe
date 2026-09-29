use qfall_math::integer::{MatZ, Z};
use qfall_math::rational::MatQ;

#[derive(Debug, Clone)]
pub struct PublicParams {
    /// computational security parameter
    pub lambda_1: usize,

    /// statistical security parameter
    pub lambda_2: usize,

    /// max rerand level
    pub t_max: usize,

    /// prime modulus q
    pub q: Z,

    /// LWE dimension n
    pub n: usize,

    /// LWE width m
    pub m: usize,

    /// k = ceil(log2(q))
    pub k: usize,

    /// LWE error distribution
    pub nu: f64,

    /// Gaussian parameter for the lattice-coset samples
    pub alpha: f64,

    /// gaussian parameter for smudging
    pub gamma: Vec<f64>,

    /// gaussian parameter for helper ciphertext smudging
    pub gamma_prime: Vec<f64>,

    /// reduced gadget basis
    pub S: MatZ,

    /// the Gram-Schmidt orthogonalization of S
    pub s_gso: MatQ,
}
