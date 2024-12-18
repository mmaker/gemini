//! The tensorcheck IOP protocol for proving multivariate evaluations
//! $f(\rho_0, \dots, \rho_{n-1})$.
//!
//! Let $f(x) \in \FF\[x\]$ be a polynomial of degree $2^n $
//! represented as the a vector of its coefficients.
//! The tensor check allows to prove the scalar product:
//!
//! \\[
//! \langle f, \otimes_j (1, \rho_j) \rangle = t
//! \\]
//!
//! for some target $t$.
//! The argument exploits even/odd folding of the polynomial.
//! That is, consider the polynomials $f_e, f_o \in \FF\[x\]$
//! of degree $2^{n-1}$:
//!
//! \\[
//! f(x) = f_e(x^2) + x f_o(x^2).
//! \\]
//!
//! Send as an oracle message
//! $ f'(x) = f_e(x) + \rho_0 f_o(x) $.
//! The verifier checks that each folded polynomial is computed correcly by
//!  testing on a random point $\beta$:
//!
//! \\[
//! f'(\beta^2) =
//!     \frac{f(\beta) + f(-\beta)}{2} + \rho_j
//!     \frac{f(\beta) - f(-\beta)}{2\beta}
//! \\]
//!
//! It proceeds recursively until the polynomial is of degree 1.
//! If we consider the map
//! $
//! \FF[x_0, \dots, x_{n-1}] \to \FF\[x\]:
//! f(x_0, \dots, x_n) \mapsto f(x, x^2, \dots, x^{2^{n-1}})
//! $
//! we are effectively
//! reducing a multivariate evaluation proof to a univariate tensorcheck.
//!
use ark_ec::pairing::Pairing;
use ark_ff::{Field, AdditiveGroup};
use ark_serialize::*;
use ark_std::borrow::Borrow;
use ark_std::vec::Vec;
use ark_std::One;

use merlin::Transcript;

use crate::iterable::Iterable;
use crate::kzg::Commitment;
use crate::kzg::CommitterKey;
use crate::kzg::EvaluationProof;
use crate::kzg::VerificationError;
use crate::kzg::VerificationResult;
use crate::kzg::VerifierKey;
use crate::misc::strip_last;
use crate::misc::{evaluate_le, fold_polynomial, ip, linear_combination, powers};
use crate::subprotocols::sumcheck::streams::FoldedPolynomialTree;
use crate::transcript::GeminiTranscript;
use crate::SPACE_TIME_THRESHOLD;

pub mod streams;

#[cfg(test)]
pub mod tests;

/// Evaluate a folded polynomial tree at the point `x`.
///
/// Make a single pass over the [`FoldedPolynomialTree`](crate::subprotocols::sumcheck::streams::FoldedPolynomialTree)
/// and return a vector storing $f^{(j)}(x)$ at the $j-1$-th position.
/// Foldings are in the interval $1, \dots, n-1$.
pub fn evaluate_folding<F, S>(polynomials: &FoldedPolynomialTree<'_, F, S>, x: F) -> Vec<F>
where
    F: Field,
    S: Iterable,
    S::Item: Borrow<F>,
{
    let mut result = vec![F::zero(); polynomials.depth()];
    for (i, c) in polynomials.iter() {
        // foldings are in the interval [[1, n-1]].
        // Align the level with the index in result.
        let i = i - 1;
        result[i] = result[i] * x + c.borrow();
    }

    result
}

/// Compute the evaluation of folded polynomials in the next round.
///
/// Let $ f(x) $ denote the polynomial in the current round,
/// $ f'(x) $ denote the folded polynomial in the next round,
/// $ \beta $ denote the evaluation point, and $ \rho $ denote the randomness for folding.
///
/// This function computes $ f'(\beta^2) = \frac{f(\beta) + f(-\beta)}{2} + \rho \cdot \frac{f(\beta) - f(-\beta)}{2\beta}$
#[inline]
pub fn evaluate_sq_fp<F: Field>(
    eval_parent_pos: &F,
    eval_parent_neg: &F,
    sp_randomness: &F,
    two_inv: &F,
    two_beta_inv: &F,
) -> F {
    (*eval_parent_pos + eval_parent_neg) * two_inv
        + (*eval_parent_pos - eval_parent_neg) * sp_randomness * two_beta_inv
}

/// The struct for the tensor check proof.
#[derive(CanonicalSerialize, PartialEq, Eq)]
pub struct TensorcheckProof<E: Pairing> {
    /// The commitments for all the folded polynomials in the tensor check.
    pub folded_polynomials_commitments: Vec<Commitment<E>>,
    /// The evaluations of all the folded polynomials in the tensor check.
    pub folded_polynomials_evaluations: Vec<[E::ScalarField; 2]>,
    /// The batched evaluation proof for both base polynomials and folded polynomials.
    pub evaluation_proof: EvaluationProof<E>,
    /// The evaluations of base polynomials, which are used to construct evaluations in the initial round of tensor check.
    pub base_polynomials_evaluations: Vec<[E::ScalarField; 3]>,
}

/// The function for folding polynomials using given challenges for each round.
/// It skips the last challenge since the result can be obtained from asserted results.
pub fn foldings_polynomial<F: Field>(polynomial: &[F], challenges: &[F]) -> Vec<Vec<F>> {
    let challenges = strip_last(challenges);
    challenges
        .iter()
        .scan(polynomial.to_vec(), |polynomial, &challenge| {
            *polynomial = fold_polynomial(polynomial, challenge);
            polynomial.clone().into()
        })
        .collect()
}

/// Store in memory all polynomial foldings after `threshold_level`
pub(crate) fn transcribe_foldings<F, S>(
    foldings: FoldedPolynomialTree<'_, F, S>,
    threshold_level: usize,
) -> Vec<Vec<F>>
where
    F: Field,
    S: Iterable,
    S::Item: Borrow<F>,
{
    let mut transcribed_foldings = vec![Vec::new(); foldings.depth() - threshold_level];
    // Filter elements in the stream that have a certain threshold and shift the index to match the vector element.
    foldings
        .iter()
        .filter_map(|(i, item)| {
            (i != 0 && i > threshold_level).then_some((i - threshold_level - 1, item))
        })
        .for_each(|(i, item)| transcribed_foldings[i].push(item));
    // Reverse each element, as they are in little endian.
    transcribed_foldings
        .iter_mut()
        .for_each(|folding| folding.reverse());
    transcribed_foldings
}

pub(crate) fn partially_foldtree<'a, F, S>(
    stream: &'a S,
    challenges: &'a [F],
) -> (FoldedPolynomialTree<'a, F, S>, Vec<Vec<F>>)
where
    F: Field,
    S: Iterable,
    S::Item: Borrow<F>,
{
    let full_foldings = FoldedPolynomialTree::new(stream, challenges);
    let threshold_level = if full_foldings.depth() > SPACE_TIME_THRESHOLD {
        full_foldings.depth() - SPACE_TIME_THRESHOLD
    } else {
        full_foldings.depth()
    };
    let transcribed_foldings = transcribe_foldings(full_foldings, threshold_level);
    let partial_foldings = FoldedPolynomialTree::new(stream, &challenges[..threshold_level]);
    (partial_foldings, transcribed_foldings)
}

impl<E: Pairing> TensorcheckProof<E> {
    /// The function for construct tensor check proof in a time-efficient way.
    ///
    /// It takes as input the randomness generator `transcript`, the committer key `ck`,
    /// the `base_polynomials` and the folded polynomials `body_polynomials`.
    ///
    /// The `base_polynomials` are part of polynomials that are evaluated at a tensor point.
    ///
    /// The folded polynomials `body_polynomials` consist of multiple tensor check intance.
    /// Each instance contains a set of folded polynomials and folding randomnesses.
    pub fn new_time<const N: usize, const M: usize>(
        transcript: &mut Transcript,
        ck: &CommitterKey<E>,
        base_polynomials: [&Vec<E::ScalarField>; N],
        body_polynomials: [(&[&Vec<E::ScalarField>], &[E::ScalarField]); M],
    ) -> TensorcheckProof<E> {
        let max_len = body_polynomials
            .iter()
            .map(|x| x.0.len())
            .fold(0, usize::max);

        let batch_challenge = transcript.get_challenge::<E::ScalarField>(b"batch_challenge");
        let batch_challenges = powers(batch_challenge, max_len);
        assert_ne!(batch_challenges.len(), 0);
        assert!(body_polynomials
            .iter()
            .all(|polynomials| polynomials.0.len() != 0));

        let batched_body_polynomials = body_polynomials.iter().map(|(polynomials, challenges)| {
            (
                linear_combination(polynomials, &batch_challenges),
                challenges,
            )
        });

        let foldings_body_polynomials = batched_body_polynomials
            .flat_map(|(polynomial, challenges)| foldings_polynomial(&polynomial, challenges))
            .collect::<Vec<_>>();
        let folded_polynomials_commitments = ck.batch_commit(&foldings_body_polynomials);

        // add commitments to transcript
        folded_polynomials_commitments
            .iter()
            .for_each(|c| transcript.append_serializable(b"commitment", c));
        let eval_chal = transcript.get_challenge::<E::ScalarField>(b"evaluation-chal");
        let minus_eval_chal = -eval_chal;
        let eval_chal2 = eval_chal.square();

        let base_polynomials_evaluations = base_polynomials
            .iter()
            .map(|polynomial| {
                [
                    evaluate_le(polynomial, &eval_chal2),
                    evaluate_le(polynomial, &eval_chal),
                    evaluate_le(polynomial, &minus_eval_chal),
                ]
            })
            .collect::<Vec<_>>();

        let folded_polynomials_evaluations = foldings_body_polynomials
            .iter()
            .map(|polynomial| {
                [
                    evaluate_le(polynomial.borrow(), &eval_chal),
                    evaluate_le(polynomial.borrow(), &minus_eval_chal),
                ]
            })
            .collect::<Vec<_>>();

        let mut all_polynomials = base_polynomials.to_vec();
        all_polynomials.extend(foldings_body_polynomials.iter());

        // add all evaluations to the transcript
        base_polynomials_evaluations
            .iter()
            .flatten()
            .for_each(|e| transcript.append_serializable(b"eval", e));
        folded_polynomials_evaluations
            .iter()
            .flatten()
            .for_each(|e| transcript.append_serializable(b"eval", e));
        let open_chal = transcript.get_challenge(b"open-chal");

        let evaluation_proof = ck.batch_open_multi_points(
            &all_polynomials[..],
            &[eval_chal2, eval_chal, minus_eval_chal],
            &open_chal,
        );

        Self {
            base_polynomials_evaluations,
            folded_polynomials_evaluations,
            evaluation_proof,
            folded_polynomials_commitments,
        }
    }

    /// The function for verifying tensor check proof.
    ///
    /// It takes as input the randomness generator `transcript`, the verifying key `vk`,
    /// the asserted result `asserted_res_vec` for each tensor check instance,
    /// the commitments `base_polynomials_commitments` for all base polynomials,
    /// the evaluations `direct_base_polynomials_evaluations` of base polynomials in each tensor check instance,
    /// folding randomnesses `fold_randomness`,
    /// the evaluation challenge point `eval_chal`,
    /// and the random challenge `batch_challenge` for batching tensor check instances with the same folding randomnesses.
    pub fn verify(
        &self,
        transcript: &mut Transcript,
        vk: &VerifierKey<E>,
        asserted_res_vec: &[Vec<E::ScalarField>],
        base_polynomials_commitments: &[Commitment<E>],
        direct_base_polynomials_evaluations: &[[E::ScalarField; 2]],
        fold_randomness: &[Vec<E::ScalarField>],
        eval_chal: E::ScalarField,
        batch_challenge: E::ScalarField,
    ) -> VerificationResult
    where
        E: Pairing,
    {
        let minus_eval_chal = -eval_chal;
        let eval_chal2 = eval_chal.square();

        let two_inv = E::ScalarField::one().double().inverse().unwrap();
        let two_beta_inv = eval_chal.double().inverse().unwrap();

        let mut evaluations = Vec::new();
        evaluations.extend(
            self.base_polynomials_evaluations
                .iter()
                .map(|x| x.to_vec())
                .collect::<Vec<_>>(),
        );

        let mut offset = 0;
        for (instance, randomness) in fold_randomness.iter().enumerate() {
            let rounds = randomness.len() - 1;
            let base_evals = &direct_base_polynomials_evaluations[instance];
            let folded_polynomials_evaluations =
                &self.folded_polynomials_evaluations[offset..offset + rounds];
            let asserted_res = &asserted_res_vec[instance];
            offset += rounds;

            evaluations.push(vec![
                evaluate_sq_fp(
                    &base_evals[0],
                    &base_evals[1],
                    &randomness[0],
                    &two_inv,
                    &two_beta_inv,
                ),
                folded_polynomials_evaluations[0][0],
                folded_polynomials_evaluations[0][1],
            ]);

            for i in 1..rounds {
                evaluations.push(vec![
                    evaluate_sq_fp(
                        &folded_polynomials_evaluations[i - 1][0],
                        &folded_polynomials_evaluations[i - 1][1],
                        &randomness[i],
                        &two_inv,
                        &two_beta_inv,
                    ),
                    folded_polynomials_evaluations[i][0],
                    folded_polynomials_evaluations[i][1],
                ]);
            }

            let subclaim = evaluate_sq_fp(
                &folded_polynomials_evaluations[rounds - 1][0],
                &folded_polynomials_evaluations[rounds - 1][1],
                &randomness[rounds],
                &two_inv,
                &two_beta_inv,
            );

            let batch_challenges = powers(batch_challenge, asserted_res.len());
            let lc_asserted_res = ip(asserted_res, &batch_challenges);

            if subclaim != lc_asserted_res {
                return Err(VerificationError);
            }
        }

        let mut all_commitments = base_polynomials_commitments.to_vec();
        all_commitments.extend(self.folded_polynomials_commitments.iter());

        self.base_polynomials_evaluations
            .iter()
            .flatten()
            .for_each(|e| transcript.append_serializable(b"eval", e));
        self.folded_polynomials_evaluations
            .iter()
            .flatten()
            .for_each(|e| transcript.append_serializable(b"eval", e));
        let open_chal = transcript.get_challenge(b"open-chal");

        vk.verify_multi_points(
            &all_commitments,
            &[eval_chal2, eval_chal, minus_eval_chal],
            &evaluations,
            &self.evaluation_proof,
            &open_chal,
        )
    }
}
#[test]
fn test_foldings_polynomial() {
    use ark_ff::One;
    use ark_test_curves::bls12_381::Fr;
    let polynomial = [Fr::from(100), Fr::from(101), Fr::from(102), Fr::from(103)];
    let challenges = [Fr::one(), Fr::one()];
    let foldings = foldings_polynomial(&polynomial, &challenges);
    assert_eq!(foldings[0].len(), 2);
    assert_eq!(foldings[0][0], Fr::from(100 + 101));
    // assert_eq!(foldings[1].len(), 1);
    // assert_eq!(foldings[1][0], Fr::from(100 + 101 + 102 + 103));
}

/// Macro rule for produce foldings for batched polynomial streams.
#[macro_export]
macro_rules! batch_polynomial_foldings {
    ($polynomials:expr, $challenges:expr) => {{
        $crate::sumcheck::streams::FoldedPolynomialTree::new(
            $crate::tensorcheck::streams::LinCombStream::new($polynomials),
            $challenges,
        )
    }};
}
