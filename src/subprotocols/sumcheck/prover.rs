//! Common data structures for the prover algorith in the scalar-product sub-argument.
use ark_ff::Field;
use ark_serialize::CanonicalSerialize;
use ark_std::boxed::Box;
use ark_std::iter::Sum;
use ark_std::vec::Vec;

/// Each message from the prover in a sumcheck protocol is a pair of FF-elements.
#[derive(CanonicalSerialize, Copy, Clone, Debug, PartialEq, Eq)]
pub struct RoundMsg<F: Field>(pub(crate) F, pub(crate) F);

/// Messages sent by the prover throughout the protocol.
#[derive(CanonicalSerialize, Clone, Debug, PartialEq, Eq)]
pub struct ProverMsgs<F: Field>(pub(crate) Vec<RoundMsg<F>>, pub(crate) Vec<[F; 2]>);

impl<F: Field> Sum for RoundMsg<F> {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.reduce(|fst, snd| RoundMsg(fst.0 + snd.0, fst.1 + snd.1))
            .unwrap_or_else(|| RoundMsg(F::zero(), F::zero()))
    }
}

impl<F: Field> RoundMsg<F> {
    pub(crate) fn mul(self, rhs: &F) -> Self {
        RoundMsg(self.0 * rhs, self.1 * rhs)
    }
}

/// Prover trait interface for both time-efficient and space-efficient prover.
pub trait Prover<F>: Send + Sync
where
    F: Field,
{
    /// Return the next prover message (if any).
    fn next_message(&mut self, verifier_message: Option<F>) -> Option<RoundMsg<F>>;
    /// Perform even/odd folding of the instance using the challenge `challenge`.
    fn fold(&mut self, challenge: F);
    // Return the total number of rouds in the protocol.
    fn rounds(&self) -> usize;
    /// Current round number.
    fn round(&self) -> usize;
    /// Return the fully-folded isntances if at the final round,
    /// otherwise return None.
    fn final_foldings(&self) -> Option<[F; 2]>;
}

impl<'a, F: Field> Prover<F> for Box<dyn Prover<F> + 'a> {
    fn next_message(&mut self, verifier_message: Option<F>) -> Option<RoundMsg<F>> {
        (**self).next_message(verifier_message)
    }

    fn fold(&mut self, challenge: F) {
        (**self).fold(challenge)
    }

    fn rounds(&self) -> usize {
        (**self).rounds()
    }

    fn round(&self) -> usize {
        (**self).round()
    }

    fn final_foldings(&self) -> Option<[F; 2]> {
        (**self).final_foldings()
    }
}
