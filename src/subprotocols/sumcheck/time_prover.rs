//!
//!

use ark_ff::Field;
use ark_std::log2;
use ark_std::vec::Vec;

use crate::misc::fold_polynomial;
use crate::subprotocols::sumcheck::prover::{Prover, RoundMsg};

/// The witness for the Twisted Scalar product relation.
#[derive(Clone)]
pub(crate) struct Witness<F: Field> {
    /// The left hand side
    f: Vec<F>,
    /// The right-hand side
    g: Vec<F>,
    /// The twist
    twist: F,
}

/// The witness for the proving algorithm.
impl<F: Field> Witness<F> {
    /// Instantiate a new Witness isntance from polynomials f, g.
    /// The instance also accepts a `twist`.
    pub(crate) fn new(f: &[F], g: &[F], twist: &F) -> Witness<F> {
        Witness {
            f: f.to_vec(),
            g: g.to_vec(),
            twist: *twist,
        }
    }

    /// Output the number of rounds required for the given scalar product.
    pub fn required_rounds(&self) -> usize {
        let max_len = usize::max(self.f.len(), self.g.len());
        log2(max_len) as usize
    }
}

/// The state of the time prover in the scalar product protocol.
pub struct TimeProver<F: Field> {
    /// The polynomial `f` in the scalar product.
    pub f: Vec<F>,
    /// The polynomial `g` in the scalar product.
    pub g: Vec<F>,
    /// Round counter.
    pub round: usize,
    pub twist: F,
    /// Total number of rounds.
    pub tot_rounds: usize,
}

impl<F: Field> TimeProver<F> {
    /// Create a new time prover.
    /// This will cause a copy of the witness.
    pub(crate) fn new(witness: Witness<F>) -> Self {
        // let twists = powers(witness.twist, witness.f.len());
        // let f = hadamard(&twists, &witness.f);
        TimeProver {
            f: witness.f.to_vec(),
            g: witness.g.to_vec(),
            round: 0usize,
            twist: witness.twist,
            tot_rounds: witness.required_rounds(),
        }
    }
}

impl<F> Prover<F> for TimeProver<F>
where
    F: Field,
{
    /// Fold the sumcheck instance (inplace).
    fn fold(&mut self, r: F) {
        // Fold the polynonomials f, g in the scalar product.
        self.f = fold_polynomial(&self.f, r * self.twist);
        self.g = fold_polynomial(&self.g, r);
        self.twist.square_in_place();
    }

    /// Time-efficient, next-message function.
    fn next_message(&mut self, verifier_message: Option<F>) -> Option<RoundMsg<F>> {
        assert!(self.round <= self.tot_rounds, "More rounds than needed.");
        // debug!("Round: {}", self.round);

        if let Some(challenge) = verifier_message {
            self.fold(challenge);
        }

        // If we already went through tot_rounds, no message must be sent.
        if self.round == self.tot_rounds {
            return None;
        }

        // Compute the polynomial of the partial sum q = a + bx + c x2,
        // For the evaluations, send only the coefficients a, b of the polynomial .
        let mut a = F::zero();
        let mut b = F::zero();
        let zero = F::zero();
        let twist2 = self.twist.square();

        let mut twist_runner = F::one();

        for (f_pair, g_pair) in self.f.chunks(2).zip(self.g.chunks(2)) {
            // The even part of the polynomial must always be unwrapped.
            let f_even = f_pair[0];
            let g_even = g_pair[0];

            // For the right part, we might obtain zero if the degree is not a multiple of 2.
            let f_odd = f_pair.get(1).unwrap_or(&zero);
            let g_odd = g_pair.get(1).unwrap_or(&zero);

            // Add to the partial sum
            a += f_even * g_even * twist_runner;
            b += (f_even * g_odd + g_even * f_odd * self.twist) * twist_runner;
            twist_runner *= twist2;
        }
        // Increment the round counter
        self.round += 1;

        Some(RoundMsg(a, b))
    }

    /// The number of rounds this prover is supposed to run on.
    #[inline]
    fn rounds(&self) -> usize {
        self.tot_rounds
    }

    fn round(&self) -> usize {
        self.round
    }

    fn final_foldings(&self) -> Option<[F; 2]> {
        (self.round == self.tot_rounds).then(|| [self.f[0], self.g[0]])
    }
}

#[test]
fn test_rounds() {
    use ark_ff::One;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_test_curves::bls12_381::Fr;
    let rng = &mut ark_std::test_rng();

    // in the smallest instance (degree-1 polynomials) a single round is necessary.
    let f = DensePolynomial::rand(1, rng);
    let g = DensePolynomial::rand(1, rng);
    let witness = Witness::<Fr>::new(&f, &g, &Fr::one());
    assert_eq!(witness.required_rounds(), 1);

    // otherwise, it should be the ceil of the log of the coefficients
    let f = DensePolynomial::rand(16, rng);
    let g = DensePolynomial::rand(1, rng);
    let witness = Witness::<Fr>::new(&f, &g, &Fr::one());
    assert_eq!(witness.required_rounds(), 5);
}

#[test]
fn test_trivial_prover() {
    use ark_ff::One;
    use ark_ff::UniformRand;
    use ark_poly::univariate::DensePolynomial;
    use ark_poly::DenseUVPolynomial;
    use ark_test_curves::bls12_381::Fr;

    let rng = &mut ark_std::test_rng();
    let one = Fr::one();

    // Instantiate a new prover for a scalar product (f . g) with f, g of degree 1.
    let f = DensePolynomial::rand(1, rng);
    let g = DensePolynomial::rand(1, rng);
    let witness = Witness::<Fr>::new(&f, &g, &one);
    let mut prover = TimeProver::new(witness);

    // Check the next-message function is correct.
    let r = Fr::rand(rng);
    let verifier_message = Some(r);
    let prover_message = prover.next_message(verifier_message);
    assert!(prover_message.is_some());

    // Check that, after forlding the instance, the witness polynomial is a constant term, and that
    // the folding operation for degree-one polynomials is identical to evaluation.
    assert_eq!(prover.f.len(), 1);
    // assert_eq!(prover.f.len()[0], prover.f[0]+ r * prover.f[1]));

    // a subsequent call to the next-message function should return None.
    assert!(prover.next_message(None).is_none());
}
