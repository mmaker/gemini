use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_std::borrow::Borrow;
use ark_std::boxed::Box;
use ark_std::vec::Vec;
use ark_std::One;

use merlin::Transcript;

use crate::iterable::Iterable;
use crate::kzg::CommitterKeyStream;
use crate::misc::evaluate_be;
use crate::subprotocols::sumcheck::{ElasticProver, Prover, SpaceProver, TimeProver};
use crate::transcript::GeminiTranscript;

use super::streams::{entry_product_streams, ProductStream, RightRotationStreamer};
use super::{EntryProduct, ProverMsgs};

impl<'a, E: Pairing, S: Iterable<Item = E::ScalarField>>
    EntryProduct<
        E,
        ElasticProver<
            SpaceProver<
                E::ScalarField,
                ProductStream<'a, E::ScalarField, S>,
                RightRotationStreamer<'a, S>,
            >,
            TimeProver<E::ScalarField>,
        >,
    >
{
    /// Create a new (single) entry product argument
    /// using the commiter key `ck` and the stream `v`, whose grand product is `claimed_product`.
    pub fn new_elastic<SG>(
        transcript: &mut Transcript,
        ck: &CommitterKeyStream<E, SG>,
        v: &'a S,
        claimed_product: E::ScalarField,
    ) -> Self
    where
        SG: Iterable,
        SG::Item: Borrow<E::G1Affine>,
    {
        let (rrot_v, acc_v) = entry_product_streams(v);

        let acc_v_commitments = vec![ck.commit(&acc_v)];
        transcript.append_serializable(b"acc_v", &acc_v_commitments[0]);

        let chal = transcript.get_challenge::<E::ScalarField>(b"ep-chal");

        let claimed_sumchecks = vec![
            chal * evaluate_be(acc_v.iter(), &chal) + claimed_product
                - chal.pow(&[acc_v.len() as u64]),
        ];
        let provers = vec![ElasticProver::new(acc_v, rrot_v, chal)];
        let msgs = ProverMsgs {
            acc_v_commitments,
            claimed_sumchecks,
        };
        EntryProduct {
            provers,
            chal,
            msgs,
        }
    }
}

macro_rules! impl_elastic_batch {
    ($name: ident; $($B:ident), *) => {
        #[allow(non_snake_case)]
        #[allow(unused_assignments)]
        pub fn $name<SG, $($B),*>(
            transcript: &mut Transcript,
            ck: &CommitterKeyStream<E, SG>,
            vs: ($(&'a $B,)*),
            claimed_products: &[E::ScalarField],
        ) -> Self
        where
            SG: Iterable,
            SG::Item: Borrow<E::G1Affine>,
            $(
                $B: crate::iterable::Iterable<Item=E::ScalarField>,
            )*

        {
            let ($($B,)*) = vs;
            let mut acc_v_commitments = Vec::new();
            $(
                let acc_v = ProductStream::new($B);
                let acc_v_commitment = ck.commit(&acc_v);
                transcript.append_serializable(b"acc_v", &acc_v_commitment);
                acc_v_commitments.push(acc_v_commitment);
            )*

            let chal = transcript.get_challenge::<E::ScalarField>(b"ep-chal");

            let mut claimed_sumchecks = Vec::new();
            let mut provers = Vec::<Box<dyn Prover<E::ScalarField> + 'a>>::new();
            let mut claimed_products_it = claimed_products.into_iter();


            $(
                let rrot_v = RightRotationStreamer::new($B, E::ScalarField::one());
                let acc_v = ProductStream::new($B);
                let claimed_product = claimed_products_it.next().expect("mismatch in claimed prod len");

                let acc_v_chal = evaluate_be(acc_v.iter(), &chal);
                let chal_n = chal.pow(&[acc_v.len() as u64]);
                let claimed_sumcheck =  acc_v_chal * chal + claimed_product - chal_n;

                claimed_sumchecks.push(claimed_sumcheck);
                let sumcheck_prover = ElasticProver::new(acc_v, rrot_v, chal);
                provers.push(Box::new(sumcheck_prover));
            )*

            let msgs = ProverMsgs {
                acc_v_commitments,
                // XXXX. should we send also the claimed sumchecks?
                claimed_sumchecks,
            };
            EntryProduct { msgs, chal, provers }
        }
    };
}

impl<'a, E: Pairing> EntryProduct<E, Box<dyn Prover<E::ScalarField> + 'a>> {
    // lets gooooo
    // impl_elastic_batch!(new_elastic_batch; A0, A1, A2);
    // impl_elastic_batch!(new_elastic_batch; A0, A1, A2, A3);
    // impl_elastic_batch!(new_elastic_batch; A0, A1, A2, A3, A4);
    // impl_elastic_batch!(new_elastic_batch; A0, A1, A2, A3, A4, A5);
    // impl_elastic_batch!(new_elastic_batch; A0, A1, A2, A3, A4, A5, A6);
    // impl_elastic_batch!(new_elastic_batch; A0, A1, A2, A3, A4, A5, A6, A7);
    impl_elastic_batch!(new_elastic_batch; A0, A1, A2, A3, A4, A5, A6, A7, A8);
}
