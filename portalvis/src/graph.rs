use crate::bits::DenseBits;
use crate::geometry::Winding;
use crate::geometry::Plane;
use std::collections::HashMap;

// following vvis convention:
//   the normal of a winding (abcdef…) is cross(a-b, c-b)
//   a prt line "src dst abdefg" has a winding that points from dst to src
//   
#[derive(Clone, Debug)]
pub struct Portal {
	pub winding: Winding,
	pub plane:   Plane,
	pub leaf_from: usize,
	pub leaf_into: usize
}

#[derive(Clone, Debug, Default)]
pub struct LeafGraph {
	pub leaf_from: Vec<Vec<usize>>, // leaf_from[leaf] = [portal, portal, portal]
	pub portals:   Vec<Portal>
}

impl LeafGraph {
	pub fn advance(&self, i: usize) -> &Vec<usize> { &self.leaf_from[self.portals[i].leaf_into] }

	pub fn retreat(&self, i: usize) -> &Vec<usize> { &self.leaf_from[self.portals[i].leaf_from] }

	pub fn filter_portals(&self, portal_indices: &Vec<usize>) -> LeafGraph {
		let mut m: HashMap<usize, usize> = HashMap::new();
		for (i, p) in portal_indices.iter().enumerate() { m.insert(*p, i); }
		LeafGraph {
			leaf_from: self.leaf_from.iter().map(|lfp| lfp.iter().flat_map(|pp| m.get(pp).copied()).collect()).collect(),
			portals:   (0..portal_indices.len()).into_iter().map(|i| { let op = &self.portals[portal_indices[i]]; Portal {
				winding: op.winding.clone(),
				plane:   op.plane,
				leaf_into: op.leaf_into,
				leaf_from: 0usize,
			}}).collect()
		}
	}

	pub fn collapse_cells(&self, ngroups: usize, cell_group: &Vec<usize>) -> LeafGraph {
		assert_eq!(cell_group.len(), self.leaf_from.len());

		let mut graph = LeafGraph {leaf_from: Vec::new(), portals: Vec::new() };
		graph.leaf_from.resize(ngroups, Vec::new());

		for Portal { winding, plane, leaf_from, leaf_into } in &self.portals {
			let group_from = cell_group[*leaf_from];
			let group_into = cell_group[*leaf_into];
			if group_from != group_into {
				graph.leaf_from[group_from].push(graph.portals.len());
				graph.portals.push(Portal {
					winding:   winding.clone(),
					plane:     *plane,
					leaf_from: group_from,
					leaf_into: group_into
				})
			}
		}

		graph
	}

	pub fn square<T, F: Fn(usize)->T>(&self, f: F) -> Vec<T> {
		let nportals = self.portals.len();
		let mut v = Vec::new();
		v.resize_with(nportals, || f(nportals));
		v
	}

	pub fn portalvis_to_leafvis(&self, portalvis: &Vec<DenseBits>) -> Vec<Vec<bool>> {
		let nportals = self.portals.len();
		let nleaves = self.leaf_from.len();
		let mut leafvis: Vec<Vec<bool>> = Vec::new();
		leafvis.resize_with(nleaves, Default::default);
		for leaf in 0..nleaves {
			leafvis[leaf].resize(nleaves, false);
		}
		for leaf in 0..nleaves {
			for p in &self.leaf_from[leaf] {
				for q in 0..nportals {
					if portalvis[*p][q] {
						leafvis[leaf][self.portals[q].leaf_into] = true;
					}
				}
				leafvis[leaf][self.portals[*p].leaf_into] = true;
			}
			leafvis[leaf][leaf] = true;
		}
		leafvis
	}
}
