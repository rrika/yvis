use std::collections::HashSet;
use crate::graph::LeafGraph;
use crate::geometry::N;

pub fn trace(graph: &LeafGraph, a: usize, b: usize, start: (N, N, N), dir: (N, N, N)) -> bool {
	let mut entered = HashSet::<usize>::new();
	let mut queue = vec![a];
	while let Some(leaf) = queue.pop() {
		if entered.contains(&leaf) { continue }
		entered.insert(leaf);
		for p in &graph.leaf_from[leaf] {
			let portal = &graph.portals[*p];
			if portal.winding.ray_passes(start, dir) {
				queue.push(portal.leaf_into)
			}
		}
	}
	entered.contains(&b)
}
