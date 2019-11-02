use crate::graph::LeafGraph;
use crate::geometry::chop_winding;
use crate::geometry::clip_hourglass;
use crate::geometry::N;
use crate::geometry::Plane;
use crate::geometry::Winding;

#[derive(Clone)]
pub struct WindingLimiter<'a>(pub Plane, pub Winding, pub Option<Winding>,
	pub Option<&'a Vec<(N, N)>>);

#[derive(Clone)]
pub struct FMLimiter<'a>(pub Vec<&'a Winding>);

#[derive(Clone)]
pub struct Unlimiter();

pub trait Limiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> where Self: Sized;
}

impl Limiter<'_> for WindingLimiter<'_> {
	fn traverse(&self, graph: &LeafGraph, portal: usize) -> Option<Self> {
		let (source_plane, mut source) = (self.0, self.1.clone());
		let (mut target_plane, mut target) = (graph.portals[portal].plane, graph.portals[portal].winding.clone());
		//println!("   target_plane = {:?}", target_plane);
		if let Some(betterdist) = self.3 {
			target_plane.3 = betterdist[portal].1;
		}
		if source_plane == -target_plane {
			//println!("   u-turn");
			return None
		}
		target = chop_winding(&target.points, source_plane)?.or(Some(target)).unwrap();
		source = chop_winding(&source.points, -target_plane)?.or(Some(source)).unwrap();
		if let Some(pass) = &self.2 {
			target = clip_hourglass(&source.points, &pass.points, &target.points)?.or(Some(target)).unwrap();
		}
		Some(WindingLimiter(source_plane, source, Some(target), self.3))
	}
}

impl<'a> Limiter<'a> for FMLimiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<FMLimiter<'a>> {
		let mut portals = self.0.clone();
		portals.push(&graph.portals[portal].winding);
		//if crate::geometry::fm_see_through_portals(&portals) {
		if crate::geometry::lpsolve_see_through_portals(&portals) {
			Some(FMLimiter(portals))
		} else {
			None
		}
	}
}

impl Limiter<'_> for Unlimiter {
	fn traverse(&self, _: &LeafGraph, _: usize) -> Option<Unlimiter> {
		Some(Unlimiter())
	}
}

impl<'a, A: Limiter<'a>, B: Limiter<'a>> Limiter<'a> for (A, B) {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> {
		let update_w = self.0.traverse(graph, portal)?;
		let update_fm = self.1.traverse(graph, portal)?;
		Some((update_w, update_fm))
	}
}
