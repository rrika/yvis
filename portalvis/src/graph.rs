use crate::geometry::Winding;
use crate::geometry::Plane;
use crate::geometry::chop_winding;
use crate::geometry::clip_hourglass;
use crate::geometry::all_below;
use crate::geometry::all_above;
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

#[derive(Debug, Default)]
pub struct LeafGraph {
	pub leaf_from: Vec<Vec<usize>>, // leaf_from[leaf] = [portal, portal, portal]
	pub leaf_into: Vec<usize>,      // leaf_into[portal] = leaf
	pub plane:     Vec<Plane>,      // plane[portal]
	pub winding:   Vec<Winding>     // winding[portal]
}

pub struct MightSeeLimiter(BitVec);
pub struct WindingLimiter(Plane, Winding, Option<Winding>);
pub struct FMLimiter<'a>(Vec<&'a Winding>);

pub trait Limiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> where Self: Sized;
}

impl<'a> Limiter<'a> for FMLimiter<'a> {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<FMLimiter<'a>> {
		let mut portals = self.0.clone();
		portals.push(&graph.winding[portal]);
		if crate::geometry::fm_see_through_portals(&portals) {
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

impl Limiter<'_> for WindingLimiter {
	fn traverse(&self, graph: &LeafGraph, portal: usize) -> Option<Self> {
		let (source_plane, mut source) = (self.0, self.1.clone());
		let (target_plane, mut target) = (graph.plane[portal], graph.winding[portal].clone());
		//println!("   target_plane = {:?}", target_plane);
		if source_plane == -target_plane {
			//println!("   u-turn");
			return None
		}
		target = chop_winding(&target, source_plane)?;
		source = chop_winding(&source, -target_plane)?;
		if let Some(pass) = &self.2 {
			target = clip_hourglass(&source, &pass, &target)?;
		}
		Some(WindingLimiter(source_plane, source, Some(target)))
	}
}

pub struct Watcher(u8, u8);
pub struct LearningLimiter(Vec<Watcher>, BitVec);


impl Limiter<'_> for LearningLimiter {
	fn traverse(&self, _graph: &LeafGraph, _portal: usize) -> Option<Self> {
		panic!()
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
			visited += recursive_leaf_flow(graph, base, graph.leaf_into[p], &new_mightsee, new_winding_limiter,
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
		flood(graph, portalfront, portalflood, p, graph.leaf_into[q]);
	}
}

pub fn process_graph(graph: &LeafGraph) {
	let nportals2 = graph.winding.len();

	println!("process_graph: init");

	let mut portalflood: Vec<BitVec> = Vec::new();
	let mut portalvis: Vec<BitVec> = Vec::new();
	let mut progress: Vec<std::sync::atomic::AtomicBool> = Vec::new();
	portalflood.resize_with(nportals2, Default::default);
	portalvis.resize_with(nportals2, Default::default);
	progress.resize_with(nportals2, || AtomicBool::new(false));

	if false {
		for i in 0..nportals2 {
			portalflood[i].resize_with(nportals2, || true);
			portalvis[i].resize_with(nportals2, Default::default);
		}
	} else {
		let mut portalfront: Vec<BitVec> = Vec::new();
		portalfront.resize_with(nportals2, Default::default);
		for i in 0..nportals2 {
			portalfront[i].resize_with(nportals2, || false);
			portalflood[i].resize_with(nportals2, || false);
			portalvis[i].resize_with(nportals2, Default::default);
		}
		println!("process_graph: portalfront");

		for i in 0..nportals2 {
			let mut c = 0;
			for j in 0..nportals2 {
				let f = i != j &&
				 	!all_below(&graph.winding[j].points, graph.plane[i]) &&
				 	!all_above(&graph.winding[i].points, graph.plane[j]);
				portalfront[i][j] = f;
				if f { c+=1; }
			}
			println!("front {:?} -> {:?}/{:?}", i, c, nportals2);
		}

		println!("process_graph: portalflood");
		for i in 0..nportals2 {
			flood(&graph, &portalfront, &mut portalflood, i, graph.leaf_into[i]);
		}
	}

	for p in 0..nportals2 {
		print!("recursive_leaf_flow {:?}/{:?}", p, nportals2);
		let mightsee = &portalflood[p];
		let limiter = WindingLimiter(
			graph.plane[p],
			graph.winding[p].clone(),
			None);
		// enable this to exhaust your RAM
		// let limiter = FMLimiter(vec![&graph.winding[p]]);
    	let visited = recursive_leaf_flow(
    		&graph,
			p,
			graph.leaf_into[p],
			&mightsee,
			limiter,
			&portalflood,
			&mut portalvis,
			&progress
    	);
		println!(" (visited {:?} leafs)", visited);
    	//println!(" vis = {:?}", portalvis[p]);
    	progress[p].store(true, Ordering::Release);
    }
}
