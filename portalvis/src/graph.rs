use crate::geometry::ndot;
use crate::geometry::Winding;
use crate::geometry::Plane;
use crate::geometry::chop_winding;
use crate::geometry::clip_hourglass;
use crate::geometry::all_below;
use crate::geometry::all_above;
use crate::geometry::N;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;


// following vvis convention:
//   the normal of a winding (abcdef…) is cross(a-b, c-b)
//   a prt line "src dst abdefg" has a winding that points from dst to src
//   

type BitVec = Vec<bool>;

pub fn bitvec_and_and_not(a: &BitVec, b: &BitVec, c: &BitVec) -> (BitVec, bool) {
	let mut o = BitVec::new();
	let mut any = false;
	o.reserve(a.len());
	for i in 0..a.len() {
		let v = a[i] && b[i] && !c[i];
		o.push(v);

		if v { any = true }
	}
	return (o, any)
}

#[derive(Debug)]
pub struct Portal {
	pub winding: Winding,
	pub plane:   Plane,
	pub leaf_from: usize,
	pub leaf_into: usize
}

#[derive(Debug, Default)]
pub struct LeafGraph {
	pub leaf_from: Vec<Vec<usize>>, // leaf_from[leaf] = [portal, portal, portal]
	pub portals:   Vec<Portal>
}

impl LeafGraph {
	pub fn advance(&self, i: usize) -> &Vec<usize> { &self.leaf_from[self.portals[i].leaf_into] }
	pub fn retreat(&self, i: usize) -> &Vec<usize> { &self.leaf_from[self.portals[i].leaf_from] }
}

pub struct MightSeeLimiter(BitVec);
#[derive(Clone)]
pub struct WindingLimiter<'a>(Plane, Winding, Option<Winding>,
	Option<&'a Vec<(N, N)>>);
#[derive(Clone)]
pub struct FMLimiter<'a>(Vec<&'a Winding>);

pub trait Limiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> where Self: Sized;
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

impl Limiter<'_> for MightSeeLimiter {
	fn traverse(&self, _graph: &LeafGraph, _portal: usize) -> Option<MightSeeLimiter> {
		panic!()
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
		target = chop_winding(&target, source_plane)?;
		source = chop_winding(&source, -target_plane)?;
		if let Some(pass) = &self.2 {
			target = clip_hourglass(&source, &pass, &target)?;
		}
		Some(WindingLimiter(source_plane, source, Some(target), self.3))
	}
}

pub struct Watcher(u8, u8);
pub struct LearningLimiter(Vec<Watcher>, BitVec);


impl Limiter<'_> for LearningLimiter {
	fn traverse(&self, _graph: &LeafGraph, _portal: usize) -> Option<Self> {
		panic!()
	}
}

pub struct CombinedLimiter<'a, 'b>(WindingLimiter<'a>, FMLimiter<'b>);

impl<'a, 'b> Limiter<'b> for CombinedLimiter<'a, 'b> {
	fn traverse(&self, graph: &'b LeafGraph, portal: usize) -> Option<Self> {
		let update_w = self.0.traverse(graph, portal)?;
		let update_fm = self.1.traverse(graph, portal)?;
		Some(CombinedLimiter(update_w, update_fm))
	}
}


pub fn recursive_leaf_flow<'a, T: Limiter<'a>>(
	graph: &'a LeafGraph,
	base: usize,
	leaf: usize,
	mightsee: &BitVec,
	winding_limiter: T,

	portalflood: &Vec<BitVec>,
	portalvis: &mut Vec<BitVec>,
	progress: &Vec<AtomicBool>) -> usize
{
	let mut visited = 0;
	for pp in &graph.leaf_from[leaf] { // in portals out of leaf
		let p = *pp;
		//println!("  {:?}", p);

		if mightsee[p] == false {
			//println!("    may not see");
			continue }

		let localsee = if progress[p].load(Ordering::Acquire) { &portalvis[p] } else { &portalflood[p] };

		let (new_mightsee, any_unexplored) = bitvec_and_and_not(
			mightsee,          // visible from base
			localsee,          // visible from p
			&portalvis[base]); // but not already confirmed

		// println!("  localsee     = {:?}", localsee);
		// println!("  new mightsee = {:?}", new_mightsee);

		if portalvis[base][p] == true && any_unexplored == false {
			//println!("    already explored");
			continue }

		if let Some(new_winding_limiter) = winding_limiter.traverse(graph, p) {
			//println!("   {{");
			portalvis[base][p] = true;
			visited += recursive_leaf_flow(graph, base, graph.portals[p].leaf_into, &new_mightsee, new_winding_limiter,
				portalflood,
				portalvis,
				progress);
			//println!("   }}");
		}
	}
	return visited+1
}

pub fn flood(
	graph: &LeafGraph,
	portalfront: &Vec<BitVec>,
	portalflood: &mut Vec<BitVec>,
	p: usize,
	leaf: usize)
{
	for rq in &graph.leaf_from[leaf] { // in portals out of leaf
		let q = *rq;
		if !portalfront[p][q] { continue }
		if portalflood[p][q] { continue }
		portalflood[p][q] = true;
		flood(graph, portalfront, portalflood, p, graph.portals[q].leaf_into);
	}
}

use std::time::Instant;

pub fn process_graph(graph: &LeafGraph) -> Vec<BitVec> {
	let nportals2 = graph.portals.len();

	println!("process_graph: init");

	let mut portalflood: Vec<BitVec> = Vec::new();
	let mut portalvis: Vec<BitVec> = Vec::new();
	let mut progress: Vec<std::sync::atomic::AtomicBool> = Vec::new();
	portalflood.resize_with(nportals2, Default::default);
	portalvis.resize_with(nportals2, Default::default);
	progress.resize_with(nportals2, || AtomicBool::new(false));

	let mut newdistances: Vec<(N, N)> = Vec::new();

	for i in 0..nportals2 {
		let plane = &graph.portals[i].plane;
		let normal = [plane.0, plane.1, plane.2];

		let ahead_dist = graph.advance(i).iter()
			.map(|j| &graph.portals[*j])
			.filter(|p| p.plane != -*plane)
			.flat_map(|p| &p.winding.points)
			.map(|pt| ndot(normal, *pt))
			.fold(std::f64::INFINITY, |a, b| a.min(b));

		let behind_dist = graph.retreat(i).iter()
			.map(|j| &graph.portals[*j])
			.filter(|p| p.plane != -*plane)
			.flat_map(|p| &p.winding.points)
			.map(|pt| ndot(normal, *pt))
			.fold(-std::f64::INFINITY, |a, b| a.max(b));

		let dd = (
			if ahead_dist == std::f64::INFINITY { plane.3 } else { ahead_dist },
			if behind_dist == -std::f64::INFINITY { plane.3 } else { behind_dist }
		);
		//println!("pplane +{:?}:-{:?}", dd.0-plane.3, plane.3-dd.1);
		newdistances.push(dd);
	}

	if false {
		for i in 0..nportals2 {
			portalflood[i].resize_with(nportals2, || true);
			portalvis[i].resize_with(nportals2, Default::default);
		}
	} else {
		let mut portalfront: Vec<BitVec> = Vec::new();
		let mut portalfront2: Vec<BitVec> = Vec::new();
		portalfront.resize_with(nportals2, Default::default);
		portalfront2.resize_with(nportals2, Default::default);
		for i in 0..nportals2 {
			portalfront[i].resize_with(nportals2, || false);
			portalfront2[i].resize_with(nportals2, || false);
			portalflood[i].resize_with(nportals2, || false);
			portalvis[i].resize_with(nportals2, Default::default);
		}
		println!("process_graph: portalfront");
		let now = Instant::now();

		for i in 0..nportals2 {
			let mut c = 0;
			let mut c2 = 0;
			for j in 0..nportals2 {
				if i == j {
					portalfront[i][j] = false;
					portalfront2[i][j] = false;
				} else {
					let f =
						!all_below(&graph.portals[j].winding.points, graph.portals[i].plane) &&
						!all_above(&graph.portals[i].winding.points, graph.portals[j].plane);
					let f2 =
						!all_below(&graph.portals[j].winding.points, graph.portals[i].plane.with_dist(newdistances[i].0)) &&
						!all_above(&graph.portals[i].winding.points, graph.portals[j].plane.with_dist(newdistances[j].1));
					if f { c+=1; }
					if f2 { c2+=1; }
					portalfront[i][j] = f;
					portalfront2[i][j] = f;
				}
			}
			/*if c2 == c {
				println!("front {:?} -> {:?}/{:?}/{:?}", i, c2, c, nportals2);
			} else {
				println!("front {:?} -> {:?}/{:?}/{:?}  -{:?}", i, c2, c, nportals2, c-c2);
			}*/
		}
		println!("process_graph: portalfront took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);


		println!("process_graph: portalflood");
		let now = Instant::now();
		for i in 0..nportals2 {
			flood(&graph, &portalfront, &mut portalflood, i, graph.portals[i].leaf_into);
		}
		println!("process_graph: portalflood took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);
		for i in 0..nportals2 {
			portalflood[i].clear();
			portalflood[i].resize_with(nportals2, || false);
		}
		let now = Instant::now();
		for i in 0..nportals2 {
			flood(&graph, &portalfront2, &mut portalflood, i, graph.portals[i].leaf_into);
		}
		println!("process_graph: portalflood2 took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);
	}

	for p in 0..nportals2 {
		print!("recursive_leaf_flow {:?}/{:?}", p, nportals2);
		let mightsee = &portalflood[p];
		let limiter_w = WindingLimiter(
			graph.portals[p].plane,
			graph.portals[p].winding.clone(),
		 	None,
		 	None);
		let limiter_w2 = WindingLimiter(
			graph.portals[p].plane.with_dist(newdistances[p].0),
			graph.portals[p].winding.clone(),
		 	None,
		 	Some(&newdistances));
		let limiter_fm = FMLimiter(vec![&graph.portals[p].winding]);
		let limiter = CombinedLimiter(limiter_w.clone(), limiter_fm.clone());
		let limiter2 = CombinedLimiter(limiter_w2.clone(), limiter_fm);
		let backup = portalvis[p].clone();
    	let visited_w = recursive_leaf_flow(
    		&graph,
			p,
			graph.portals[p].leaf_into,
			&mightsee,
			limiter_w,
			&portalflood,
			&mut portalvis,
			&progress
    	);
    	portalvis[p] = backup.clone();
    	let visited_w2 = recursive_leaf_flow(
    		&graph,
			p,
			graph.portals[p].leaf_into,
			&mightsee,
			limiter_w2,
			&portalflood,
			&mut portalvis,
			&progress
    	);
    	portalvis[p] = backup.clone();
    	let visited = recursive_leaf_flow(
    		&graph,
			p,
			graph.portals[p].leaf_into,
			&mightsee,
			limiter,
			&portalflood,
			&mut portalvis,
			&progress
    	);
    	portalvis[p] = backup;
    	let visited2 = recursive_leaf_flow(
    		&graph,
			p,
			graph.portals[p].leaf_into,
			&mightsee,
			limiter2,
			&portalflood,
			&mut portalvis,
			&progress
    	);
    	// compare: winding
    	// against: winding with improvements
    	// against: winding + linear programming
    	// against: winding with improvements + linear programming (should be same as previous)
    	if visited_w2 != visited_w {
    		print!("   ");
    	}
		println!(" (visited {:?}/{:?}/{:?}/{:?} leafs)", visited2, visited, visited_w2, visited_w);
		// println!(" (visited {:?} leafs)", visited_w);
    	// println!(" vis = {:?}", portalvis[p]);
    	progress[p].store(true, Ordering::Release);
    }


    // todo: downsample from leaf to portal vis to leaf to leaf vis
    portalvis
}
