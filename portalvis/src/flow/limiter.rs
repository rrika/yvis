use std::fmt::Debug;
use crate::geometry::generalized_clip_hourglass;
use crate::geometry::N;
use crate::geometry::Plane;
use crate::geometry::PlaneChoppable;
use crate::geometry::PlaneSeparable;
use crate::graph::LeafGraph;
use crate::winding::chop_winding;
use crate::winding::clip_hourglass;
use crate::winding::lpsolve_see_through_portals;
use crate::winding::Winding;

#[derive(Clone, Debug)]
pub struct WindingLimiter<'a>(pub Plane, pub Winding, pub Option<Winding>,
	pub Option<&'a Vec<(N, N)>>);

#[derive(Clone, Debug)]
pub struct ChopLimiter<T: PlaneSeparable<T> + PlaneChoppable>(
	pub Plane,
	pub T,
	pub Option<T>);

#[derive(Clone)]
pub struct FMLimiter<'a>(pub Vec<&'a Winding>);

#[derive(Clone)]
pub struct Unlimiter();

pub trait Limiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> where Self: Sized;
	fn stats(&self) {
		println!("");
	}
}

impl <T:
	PlaneSeparable<T> +
	PlaneChoppable<Output=T> +
	for<'a > From<&'a Winding> +
	Debug +
	Clone
> Limiter<'_> for ChopLimiter<T> {
	fn traverse(&self, graph: &LeafGraph, portal: usize) -> Option<Self> {
		let (source_plane, mut source) = (self.0, self.1.clone());
		let target_portal = &graph.portals[portal];
		let target_plane = target_portal.plane;
		if source_plane == -target_plane {
			//println!("   u-turn");
			return None
		}

		let mut target: T = (&target_portal.winding).into();
		target = target.chop(source_plane)?.or(Some(target)).unwrap();
		source = source.chop(-target_plane)?.or(Some(source)).unwrap();
		if let Some(pass) = &self.2 {
			target = generalized_clip_hourglass(&source, &pass, &target)?.or(Some(target)).unwrap();
		}
		Some(ChopLimiter(source_plane, source, Some(target)))
	}
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
	fn stats(&self) {
		if let Some(two) = &self.2 {
			println!("a={} b={}", self.1.points.len(), two.points.len());
		} else {
			println!("a={}", self.1.points.len());
		}
	}
}

impl<'a> Limiter<'a> for FMLimiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<FMLimiter<'a>> {
		let mut portals = self.0.clone();
		portals.push(&graph.portals[portal].winding);
		//if crate::geometry::fm_see_through_portals(&portals) {
		if lpsolve_see_through_portals(&portals) {
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
		let update_a = self.0.traverse(graph, portal)?;
		let update_b = self.1.traverse(graph, portal)?;
		Some((update_a, update_b))
	}
}

#[derive(Copy, Clone)]
pub struct CheckLimiter<A, B>(pub A, pub Option<B>);

impl<'a, A: Limiter<'a> + Debug, B: Limiter<'a> + Debug> Limiter<'a> for CheckLimiter<A, B> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> {
		let update_a = self.0.traverse(graph, portal);
		if let Some(self1) = &self.1 {
			let update_b = self1.traverse(graph, portal);
			match (update_a, update_b) {
				(None,    None   ) => None,
				(None,    Some(_)) => {
					println!("{:?}", graph.portals[portal]);
					println!("{:?}", self.0);
					println!("{:?}", self.1);
					panic!("rejected by limiter A")
				},
				(Some(a), None   ) => Some(CheckLimiter(a, None)),
				(Some(a), Some(b)) => Some(CheckLimiter(a, Some(b)))
			}
		} else {
			Some(CheckLimiter(update_a?, None))
		}
	}
}