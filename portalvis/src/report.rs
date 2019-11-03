
use crate::flow::matrix_make_symmetric;
use crate::flow::calc_portalfront;
use crate::flow::calc_portalflood;
use crate::flow::calc_portalvis;
use crate::graph::LeafGraph;
use std::time::Instant;
use crate::bits::DenseBits;

#[derive(Debug)]
pub struct Traversal {
	pub approx: usize,
	pub steps:  usize,
	pub result: usize,
	pub wasted: usize,
	pub depth:  usize
}

pub fn report<'a>(graph: &'a LeafGraph) -> (Vec<DenseBits>, Vec<DenseBits>) {
	let (portalflood, ref_portalvis) = {
		let mkbitvec = |n| { let mut v = Vec::new(); v.resize((n+63)/64, 0u64); DenseBits(v) };
		let mut portalfront: Vec<DenseBits> = graph.square(mkbitvec);
		let mut portalflood: Vec<DenseBits> = graph.square(mkbitvec);
		let mut portalvis: Vec<DenseBits> = graph.square(mkbitvec);
		let now = Instant::now();
		calc_portalfront(&graph, &mut portalfront, None);
		calc_portalflood(&graph, &mut portalflood, &portalfront);
		matrix_make_symmetric(&mut portalflood);
		calc_portalvis(&graph, &mut portalvis, &portalflood);
		matrix_make_symmetric(&mut portalvis);
		println!("report: baseline took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);
		(portalflood, portalvis)
	};
	(portalflood, ref_portalvis)
}
