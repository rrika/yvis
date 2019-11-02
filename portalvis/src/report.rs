
use crate::graph::LeafGraph;

#[derive(Debug)]
pub struct Traversal {
	pub approx: usize,
	pub steps:  usize,
	pub result: usize,
	pub depth:  usize
}

pub fn report<'a>(_: &'a LeafGraph) {

}
