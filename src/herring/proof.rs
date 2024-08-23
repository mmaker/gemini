// //! Scalar-product proof implementation.
use ark_std::boxed::Box;
use ark_std::vec::Vec;
use ark_std::Zero;
use merlin::Transcript;

#[cfg(feature = "parallel")]
use rayon::iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator};

use crate::herring::prover::SumcheckMsg;
use crate::herring::Prover;
use crate::transcript::GeminiTranscript;

use super::module::BilinearModule;
use super::subclaim::SumcheckMsgs;
use super::time_prover::Witness;
use super::TimeProver;

/// A scalar product proof, containing non-oracle messages, and oracle messages together with their queries and evaluations.
#[derive(Debug, PartialEq, Eq)]
pub struct Sumcheck<M: BilinearModule> {
    /// The non-oracle messages sent throughout the protocol.
    pub messages: Vec<SumcheckMsg<M::Target>>,
    /// The challenges sent throughout the protocol.
    pub challenges: Vec<M::ScalarField>,
    /// The number of rounds in the protocol.
    pub(crate) rounds: usize,
    // Folded statements
    pub(crate) final_foldings: Vec<(M::Lhs, M::Rhs)>,
}

impl<M: BilinearModule> Sumcheck<M> {
    /// Prove function for the scalar product.
    /// The input contains a randomness generator and a prover struct.
    /// The prover struct can be either time-efficient or space-efficient
    /// depending on the configuration.
    pub fn prove<P: Prover<M>>(transcript: &mut Transcript, mut prover: P) -> Self {
        let rounds = prover.rounds();
        let mut messages = Vec::with_capacity(rounds);
        let mut challenges = Vec::with_capacity(rounds);

        let mut verifier_message = None;
        while let Some(message) = prover.next_message(verifier_message) {
            // add the message sent to the transcript
            transcript.append_serializable(b"evaluations", &message);
            // compute the challenge for the next round
            let challenge = transcript.get_challenge::<M::ScalarField>(b"challenge");
            verifier_message = Some(challenge);

            // add the message to the final proof
            messages.push(message);
            challenges.push(challenge);
        }

        let rounds = prover.rounds();
        let final_foldings = vec![prover.final_foldings().unwrap()];
        // Add the final foldings to the transcript
        transcript.append_serializable(b"final-folding", &final_foldings[0].0);
        transcript.append_serializable(b"final-folding", &final_foldings[0].1);

        Sumcheck {
            messages,
            challenges,
            rounds,
            final_foldings,
        }
    }

    /// Prove function for a batch of scalar product instances.
    pub fn prove_batch<'a>(
        transcript: &mut Transcript,
        mut provers: Vec<Box<dyn Prover<M> + 'a>>,
    ) -> Sumcheck<M> {
        let rounds = provers.iter().map(|p| p.rounds()).fold(0, usize::max);
        let mut messages: Vec<SumcheckMsg<M::Target>> = Vec::with_capacity(rounds);
        let mut challenges = Vec::with_capacity(rounds);

        let coefficients = (0..provers.len())
            .map(|_| transcript.get_challenge::<M::ScalarField>(b"batch-sumcheck"))
            .collect::<Vec<_>>();

        let mut verifier_message = None;
        for _ in 0..rounds {
            // obtain the next message from each prover, if possible in parallel.
            let round_messages = cfg_iter_mut!(provers).map(|p| {
                p.next_message(verifier_message).unwrap_or_else(|| {
                    let final_foldings = p.final_foldings().expect(
                        "If next_message is None, we expect final foldings to be available",
                    );
                    SumcheckMsg(M::p(final_foldings.0, final_foldings.1), M::Target::zero())
                })
            });
            // compute the non-oracle messagein the sumcheck:
            let message = round_messages
                .zip(&coefficients) // take the combination of messages and coefficients
                .map(|(m, c)| m * c) // multiply them if there's an actual message
                .sum(); // finally, add them up.
            messages.push(message); // add the message sent to the transcript
            transcript.append_serializable(b"evaluations", &message);
            // compute the challenge for the next round
            let challenge = transcript.get_challenge(b"challenge");
            verifier_message = Some(challenge);
            challenges.push(challenge);
        }
        let final_foldings = provers
            .iter()
            .map(|p| {
                let final_folding = p.final_foldings().unwrap();
                transcript.append_serializable(b"final-folding-lhs", &final_folding.0);
                transcript.append_serializable(b"final-folding-rhs", &final_folding.1);
                final_folding
            })
            .collect::<Vec<_>>();
        Sumcheck {
            messages,
            challenges,
            rounds,
            final_foldings,
        }
    }

    /// Create a new Proof using the Time prover.
    pub fn new_time(
        transcript: &mut Transcript,
        f: &[M::Lhs],
        g: &[M::Rhs],
        twist: &M::ScalarField,
    ) -> Self {
        let witness = Witness::new(f, g, twist);
        let prover = TimeProver::new(witness);

        Self::prove(transcript, prover)
    }

    /// Return the prover's messages.
    pub fn prover_messages(&self) -> SumcheckMsgs<M> {
        SumcheckMsgs(self.messages.clone(), self.final_foldings.clone())
    }
}
