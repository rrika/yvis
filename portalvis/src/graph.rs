use crate::geometry::ndot;
use crate::geometry::Winding;
use crate::geometry::Plane;
use crate::geometry::chop_winding;
use crate::geometry::clip_hourglass;
use crate::geometry::all_below;
use crate::geometry::all_above;
use crate::geometry::N;
use std::sync::atomic::AtomicU64;
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

pub struct CombinedLimiter<'a, 'b>(WindingLimiter<'a>, FMLimiter<'b>);

impl<'a, 'b> Limiter<'b> for CombinedLimiter<'a, 'b> {
	fn traverse(&self, graph: &'b LeafGraph, portal: usize) -> Option<Self> {
		let update_w = self.0.traverse(graph, portal)?;
		let update_fm = self.1.traverse(graph, portal)?;
		Some(CombinedLimiter(update_w, update_fm))
	}
}
impl<'a, A: Limiter<'a>, B: Limiter<'a>> Limiter<'a> for (A, B) {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> {
		let update_w = self.0.traverse(graph, portal)?;
		let update_fm = self.1.traverse(graph, portal)?;
		Some((update_w, update_fm))
	}
}


pub fn recursive_leaf_flow<'a, T: Limiter<'a>>(
	graph: &'a LeafGraph,
	work: &Work,
	winding_limiter: T,
	base: usize,

	confirmsee: &mut BitVec,
	leaf: usize,
	mightsee: &BitVec
) -> usize {
	let mut visited = 0;
	for pp in &graph.leaf_from[leaf] { // in portals out of leaf
		let p = *pp;
		//println!("  {:?}", p);

		if mightsee[p] == false {
			//println!("    may not see");
			continue }

		let (new_mightsee, any_unexplored) = {
			let guard = &epoch::pin();
			let localsee = work.get_row(p, guard);
			bitvec_and_and_not(
				mightsee,          // visible from base
				localsee,          // visible from p
				&confirmsee) // but not already confirmed
		};

		// println!("  localsee     = {:?}", localsee);
		// println!("  new mightsee = {:?}", new_mightsee);

		if confirmsee[p] == true && any_unexplored == false {
			//println!("    already explored");
			continue }

		if let Some(new_winding_limiter) = winding_limiter.traverse(graph, p) {
			//println!("   {{");
			confirmsee[p] = true;
			visited += recursive_leaf_flow(
				graph,
				&work,
				new_winding_limiter,
				base,
				confirmsee,
				graph.portals[p].leaf_into,
				&new_mightsee,
			);
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


pub fn recursive_leaf_flow_simple<'a, T: Limiter<'a>>(
	graph: &'a LeafGraph,
	work: &Work,
	limiter: T,
	base: usize
) -> (BitVec, usize, usize)
{
	let mut basevis: BitVec = Vec::new();
	basevis.resize_with(graph.portals.len(), Default::default);
	let mightsee: BitVec = {
		let guard = &epoch::pin();
		work.get_row(base, guard).clone()
	};
	let nvisited = recursive_leaf_flow(
		&graph,
		&work,
		limiter,
		base,
		&mut basevis,
		graph.portals[base].leaf_into,
		&mightsee
    );
    let mut nvisible = 0;
    for target in &basevis {
    	if *target { nvisible += 1 }
    }
    (basevis, nvisible, nvisited-1)
}

use std::time::Instant;
use std::sync::Arc;

use crossbeam::thread;
use crossbeam::epoch::{self as epoch, Atomic, Owned, Shared, Guard};
use std::borrow::Borrow;

pub struct Work {
	approx: Vec<Atomic<BitVec>>,
	first_time:     AtomicU64,
	second_time:    AtomicU64,
	first_visible:  AtomicU64,
	second_visible: AtomicU64,
	first_visited:  AtomicU64,
	second_visited: AtomicU64,
}

impl Work {
	fn get_row<'g>(&self, a: usize, guard: &'g Guard) -> &'g BitVec {
		unsafe { self.approx[a].load(Ordering::SeqCst, guard).as_ref() }.expect("why is any row null?")
	}
	fn set_row(&self, a: usize, data: BitVec) {
		self.approx[a].store(Owned::new(data), Ordering::SeqCst)
	}
}

pub fn process_graph<'a>(zgraph: &'a LeafGraph) -> (Vec<BitVec>, Vec<BitVec>) {
	let nportals2 = zgraph.portals.len();

	println!("process_graph: init");

	let mut portalflood: Vec<BitVec> = Vec::new();
	portalflood.resize_with(nportals2, Default::default);

	let mut newdistances: Vec<(N, N)> = Vec::new();

	for i in 0..nportals2 {
		let plane = &zgraph.portals[i].plane;
		let normal = [plane.0, plane.1, plane.2];

		let ahead_dist = zgraph.advance(i).iter()
			.map(|j| &zgraph.portals[*j])
			.filter(|p| p.plane != -*plane)
			.flat_map(|p| &p.winding.points)
			.map(|pt| ndot(normal, *pt))
			.fold(std::f64::INFINITY, |a, b| a.min(b));

		let behind_dist = zgraph.retreat(i).iter()
			.map(|j| &zgraph.portals[*j])
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
						!all_below(&zgraph.portals[j].winding.points, zgraph.portals[i].plane) &&
						!all_above(&zgraph.portals[i].winding.points, zgraph.portals[j].plane);
					let f2 =
						!all_below(&zgraph.portals[j].winding.points, zgraph.portals[i].plane.with_dist(newdistances[i].0)) &&
						!all_above(&zgraph.portals[i].winding.points, zgraph.portals[j].plane.with_dist(newdistances[j].1));
					if f { c+=1; }
					if f2 { c2+=1; }
					portalfront[i][j] = f;
					portalfront2[i][j] = f;
				}
			}
			if c2 == c {
				//println!("front {:?} -> {:?}/{:?}/{:?}", i, c2, c, nportals2);
			} else {
				//println!("front {:?} -> {:?}/{:?}/{:?}  -{:?}", i, c2, c, nportals2, c-c2);
			}
		}
		println!("process_graph: portalfront took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);


		println!("process_graph: portalflood");
		let now = Instant::now();
		for i in 0..nportals2 {
			flood(&zgraph, &portalfront, &mut portalflood, i, zgraph.portals[i].leaf_into);
		}
		println!("process_graph: portalflood took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);
		for i in 0..nportals2 {
			portalflood[i].clear();
			portalflood[i].resize_with(nportals2, || false);
		}
		let now = Instant::now();
		for i in 0..nportals2 {
			flood(&zgraph, &portalfront2, &mut portalflood, i, zgraph.portals[i].leaf_into);
		}
		println!("process_graph: portalflood2 took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);
	}

	let portalvis: Vec<Atomic<BitVec>> =
		portalflood.iter().map(|row| Atomic::new(row.clone())).collect();

	let work = Work {
		approx: portalvis,
		first_time: 0.into(),
		second_time: 0.into(),
		first_visible: 0.into(),
		second_visible: 0.into(),
		first_visited: 0.into(),
		second_visited: 0.into(),
	};
	let graph: Arc<_> = Arc::new(zgraph.clone());
	let graph2 = graph.clone();
	let portalflood2 = portalflood.clone();

	let workunit = move |work: &Work, p: usize| {
		print!("recursive_leaf_flow {:?}/{:?}", p, nportals2);
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
		let _limiter = CombinedLimiter(limiter_w.clone(), limiter_fm.clone());
		let _limiter2 = CombinedLimiter(limiter_w2.clone(), limiter_fm);

		let mut cc = 0;
		// for v in work.get_row(p, guard).load() {
		// 	if *v { cc+=1 }
		// }

		let now = Instant::now();
    	let (row, visible_w, visited_w) = recursive_leaf_flow_simple(
    		&graph, &work, limiter_w.clone(), p);
    	let first_pass = now.elapsed().as_micros();
    	// let (_, visited_w2) = recursive_leaf_flow_simple(
    	// 	&graph, &woget_rowrk, limiter_w2, p);
    	// let (row, visited) = recursive_leaf_flow_simple(
    	// 	&graph, &work, limiter, p);
    	// let (_, visited2) = recursive_leaf_flow_simple(
    	// 	&graph, &work, limiter2, p);

    	work.set_row(p, row);

		let now = Instant::now();
    	let (_, second_visible_w, second_visited_w) = recursive_leaf_flow_simple(
    		&graph, &work, limiter_w, p);
    	let second_pass = now.elapsed().as_micros();

    	work.first_time.fetch_add(first_pass as u64, Ordering::SeqCst);
    	work.first_visible.fetch_add(visible_w as u64, Ordering::SeqCst);
    	work.first_visited.fetch_add(visited_w as u64, Ordering::SeqCst);
    	work.second_time.fetch_add(second_pass as u64, Ordering::SeqCst);
    	work.second_visible.fetch_add(second_visible_w as u64, Ordering::SeqCst);
    	work.second_visited.fetch_add(second_visited_w as u64, Ordering::SeqCst);

    	// compare: winding
    	// against: winding with improvements
    	// against: winding + linear programming
    	// against: winding with improvements + linear programming (should be same as previous)
		//println!(" (visited {:?}/{:?}/{:?}/{:?} leafs)", visited2, visited, visited_w2, visited_w);
		if visited_w != visible_w {
			println!(" (visited {:?}/{:?} leafs, of which {:?} are visible)", visited_w, cc, visible_w);
		} else {
			println!(" (visited {:?}/{:?} leafs)", visited_w, cc);
		}
    	// println!(" vis = {:?}", portalvis[p]);
    };
	// for p in 0..nportals2 {
	// 	workunit(p);
	// }
	thread::scope(|s| {
		let rwork = &work;
		let rworkunit = &workunit;
		for i in 0..4 {
			s.spawn(move |_| {
			 	for p in 0..nportals2 {
			 		if p % 4 == i {
						rworkunit(rwork, p);
					}
				}
			});
		}
	}).unwrap();

	println!("{} {} {} {} {}",
		work.first_time.load(Ordering::Acquire),
		work.second_time.load(Ordering::Acquire),
		work.first_visible.load(Ordering::Acquire),
		work.first_visited.load(Ordering::Acquire) - work.first_visible.load(Ordering::Acquire),
		work.second_visited.load(Ordering::Acquire) - work.second_visible.load(Ordering::Acquire)
	);

	let owned_rows: Vec<Box<BitVec>> = work.approx.into_iter().map(|atomic|
		unsafe { atomic.load(Ordering::Acquire, &epoch::pin()).into_owned().into_box() }).collect();

	let portalvis: Vec<BitVec> = owned_rows.into_iter().map(|owned| *owned).collect();

	(portalvis_to_leafvis(&graph2, &portalflood2), portalvis_to_leafvis(&graph2, &portalvis))
}

fn portalvis_to_leafvis(graph: &LeafGraph, portalvis: &Vec<BitVec>) -> Vec<BitVec> {
    let nleaves = graph.leaf_from.len();
	let mut leafvis: Vec<BitVec> = Vec::new();
	leafvis.resize_with(nleaves, Default::default);
	for leaf in 0..nleaves {
		leafvis[leaf].resize(nleaves, false);
	}
	for leaf in 0..nleaves {
		for p in &graph.leaf_from[leaf] {
			for q in 0..portalvis[*p].len() {
				if portalvis[*p][q] {
					leafvis[leaf][graph.portals[q].leaf_into] = true;
				}
			}
		}
		leafvis[leaf][leaf] = true;
	}
    leafvis
}
