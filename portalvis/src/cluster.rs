pub trait ClusterProcess {
	type Item: Clone;
	type ItemEx: Clone + Ord;
	type Cluster: Clone; // also, must support *a|*b which isn't expressible in current rust
	type ClusterEx: Clone;

	fn unit(i: &Self::Item) -> Self::Cluster;

	fn ex(i: &Self::Item) -> Self::ItemEx;
	fn ex_score<'a>(ix: &'a Self::ItemEx) -> usize;
	fn cl<'a>(ix: &'a Self::ItemEx) -> Self::ClusterEx;
	fn cl_score<'a>(ex: &'a Self::ClusterEx) -> usize;

	fn add_ex(c: &Self::Cluster, cx: &Self::ClusterEx, i: &Self::Item, ix: &Self::ItemEx) -> Self::ClusterEx;
	fn add(c: &mut Self::Cluster, i: &Self::Item, ix: &Self::ItemEx);
}

use crate::bits::bitvec_and_not_pop;
use crate::bits::DenseBits;

pub struct DenseVisClusterProcess();

impl ClusterProcess for DenseVisClusterProcess {

	type Item = DenseBits;
	type Cluster = DenseBits;
	type ItemEx = usize;
	type ClusterEx = (usize, usize, usize);

	fn unit(i: &DenseBits) -> DenseBits {
		i.clone()
	}

	fn ex(i: &DenseBits) -> usize {
		i.count_ones()
	}
	fn ex_score(ix: &usize) -> usize {
		Self::cl_score(&Self::cl(ix))
	}

	fn cl(ix: &usize) -> Self::ClusterEx {
		(*ix, *ix, 1)
	}
	fn cl_score(cx: &Self::ClusterEx) -> usize {
		let (union_pop, total_pop, nmembers) = *cx;
		let precompute_steps = if nmembers > 1 {
			nmembers * union_pop
		} else { 0 };
		let main_steps = total_pop;
		let steps = precompute_steps + main_steps;
		let step_matrix_work = union_pop * union_pop;
		let work = steps * step_matrix_work;
		work / nmembers
	}
	fn add_ex(
		c:  &DenseBits,
		cx: &Self::ClusterEx,
		i:  &DenseBits,
		ix: &usize) -> Self::ClusterEx
	{(
		cx.0 + bitvec_and_not_pop(c, i),
		cx.1 + *ix,
		cx.2 + 1
	)}
	fn add(
		c:  &mut DenseBits,
		i:  &DenseBits,
		_ix: &usize)
	{
		*c = &*c | i;
	}
}

use std::collections::BinaryHeap;

pub fn greedy<T: ClusterProcess>(items: &[&T::Item]) -> Vec<Vec<usize>> {

	let tiebreaker = 1;

	let mut cchoices = Vec::new();
	let mut h: BinaryHeap<(usize, usize, T::ItemEx)> = BinaryHeap::new();
	for i in 0..items.len() {
		let ex = T::ex(&items[i]);
		let ex_score = T::ex_score(&ex);
		//if ip > 0 {
		h.push((ex_score, i, ex))
		//}
	}

	while h.len() > 0 {
		let (_, initial, initialx) = h.pop().unwrap();
		let mut cluster_indices: Vec<usize> = Vec::new();
		let mut cluster_content = T::unit(&items[initial]);
		let mut cluster_x: T::ClusterEx = T::cl(&initialx);
		let mut cluster_score = T::cl_score(&cluster_x);

		// oc = old candidate score
		// nc = new candidate score
		while let Some((oc, index, x)) = h.pop() {
			//let newuniques = bitvec_and_not_pop(items[index], &cbv);
			let ncx = T::add_ex(&cluster_content, &cluster_x, items[index], &x);
			let nc = T::cl_score(&ncx);
			// not sound when underestimating items
			assert!(nc <= oc);

			if nc < oc {
				// disappointed by this candidate
				// send them back to the queue
				h.push((nc, index, x));
				continue

			} else if  nc > oc {
				// underestimated this candidate
				panic!("rating function not sound");

			} else if nc == oc {
				// the best the queue has to offer
				// but is it enough?
				if nc >= cluster_score + tiebreaker {
					cluster_indices.push(index);
					T::add(&mut cluster_content, items[index], &x);
					cluster_x = ncx;
					cluster_score = nc;
				} else {
					break;
				}
			}
		}
		cchoices.push(cluster_indices);

		// reset all scores
		let mut h2 = BinaryHeap::new();
		for (_, i, ex) in h {
			let ex_score = T::ex_score(&ex);
			h2.push((ex_score, i, ex));
		}
		h = h2;
	}

	cchoices
}
