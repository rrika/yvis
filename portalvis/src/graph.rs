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

use std::collections::HashMap;
use std::collections::HashSet;

type BitVec = Vec<bool>;


// following vvis convention:
//   the normal of a winding (abcdef…) is cross(a-b, c-b)
//   a prt line "src dst abdefg" has a winding that points from dst to src
//   

pub fn bitvec_pop(a: &BitVec) -> usize {
	let mut c = 0;
	for i in 0..a.len() {
		if a[i] { c+=1 }
	}
	c
}
pub fn bitvec_and_not_pop(a: &BitVec, b: &BitVec) -> usize {
	let mut c = 0;
	for i in 0..a.len() {
		if a[i] && !b[i] { c+=1 }
	}
	c
}

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

pub fn bitvec_or(a: &BitVec, b: &BitVec) -> BitVec {
	let mut o = BitVec::new();
	o.reserve(a.len());
	for i in 0..a.len() {
		let v = a[i] || b[i];
		o.push(v);
	}
	o
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

fn renumber(g: &LeafGraph, rv: &Vec<usize>) -> LeafGraph {
	let mut m: HashMap<usize, usize> = HashMap::new();
	for (i, p) in rv.iter().enumerate() { m.insert(*p, i); }
	LeafGraph {
		leaf_from: g.leaf_from.iter().map(|lfp| lfp.iter().flat_map(|pp| m.get(pp).copied()).collect()).collect(),
		portals:   (0..rv.len()).into_iter().map(|i| { let op = &g.portals[rv[i]]; Portal {
			winding: op.winding.clone(),
			plane:   op.plane,
			leaf_into: op.leaf_into,
			leaf_from: 0usize,
		}}).collect()
	}
}

fn reduce_matrix(m: &Vec<BitVec>, r: &Vec<usize>) -> Vec<BitVec> {
	r.iter().map(|x|
		r.iter().map(|y|
			m[*x][*y]
		).collect::<BitVec>()
	).collect()
}

fn reinject_matrix(m: &mut Vec<BitVec>, n: Vec<BitVec>, p: &Vec<usize>) {
	for i in 0..p.len() {
		for j in 0..p.len() {
			m[p[i]][p[j]] = n[i][j]
		}
	}
}

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

impl<'a, A: Limiter<'a>, B: Limiter<'a>> Limiter<'a> for (A, B) {
	fn traverse(&self, graph: &'a LeafGraph, portal: usize) -> Option<Self> {
		let update_w = self.0.traverse(graph, portal)?;
		let update_fm = self.1.traverse(graph, portal)?;
		Some((update_w, update_fm))
	}
}

use std::borrow::Cow;
use std::ops::{Index};

pub trait BitVecLike: Sized + Index<usize, Output=bool> {}
impl<T> BitVecLike for T where T: Sized + Index<usize, Output=bool> {}

#[derive(Clone)]
pub enum MaybeSparse {
	Sparse(Vec<usize>),
	Dense(BitVec)
}

// always   3566.888ms
// >= 400	4157.697ms
// >= 300   4380.258ms
// >= 200	4712.477ms
// >= 150   5074.123ms
// >= 100	5468.178ms
// >=  50   6652.418ms

pub fn maybe_convert(bitvec: &BitVec) -> Option<Vec<usize>> {
	let mut v = Vec::<usize>::new();
	for i in 0..bitvec.len() {
		if bitvec[i] {
			// if v.len() >= 400 {
			// 	return None
			// }
			v.push(i)
		}
	}
	Some(v)
}

const SOMETRUE: bool = true;
const SOMEFALSE: bool = false;

impl Index<usize> for MaybeSparse {
	type Output = bool;
	fn index(&self, index: usize) -> &bool {
		match self {
			MaybeSparse::Sparse(values) => if values.binary_search(&index).is_ok() { &SOMETRUE } else { &SOMEFALSE },
			MaybeSparse::Dense(bitvec) => &bitvec[index]
		}
	}
}

pub trait Mightsee<Filter, Confirm> : BitVecLike + Sized + Clone where
	Filter: BitVecLike,
	Confirm: BitVecLike
{
	fn filter_for_unexplored<'a>(&'a self, f: &Filter, c: &Confirm) -> (Cow<'a, Self>, bool);
}

impl Mightsee<BitVec, BitVec> for BitVec {
	fn filter_for_unexplored<'a>(&'a self, f: &BitVec, c: &BitVec) -> (Cow<'a, BitVec>, bool) {
		//let p_self = bitvec_pop(self);
		//let p_f = bitvec_pop(f);
		//let p_c = bitvec_pop(c);
		//println!("&&~ {} {} {}", p_self, p_f, p_c);
		let (nmc, any_unex) = bitvec_and_and_not(self, f, c);
		(Cow::Owned(nmc), any_unex)
	}
}

impl Mightsee<BitVec, BitVec> for MaybeSparse {
	fn filter_for_unexplored<'a>(&'a self, f: &BitVec, c: &BitVec) -> (Cow<'a, MaybeSparse>, bool) {
		use MaybeSparse::*;
		match self {
			Sparse(values) => {
				//println!("&&~ sparse {}", values.len());
				let newvalues: Vec<usize>
					= values.iter().map(|v|*v).filter(|i: &usize| f[*i] && !c[*i]).collect();
				let any_unexplored = newvalues.len() > 0;
				(Cow::Owned(Sparse(newvalues)), any_unexplored)
			},
			Dense(bitvec) => {
				let (r, any_unex) = bitvec.filter_for_unexplored(f, c);
				match maybe_convert(&r) {
					Some(newvalues) => (Cow::Owned(Sparse(newvalues)), any_unex),
					None            => (Cow::Owned(Dense(r.into_owned())), any_unex)
				}
			}
		}
	}
}

pub fn recursive_leaf_flow<
	'a,
	T: Limiter<'a>,
	Filter:  BitVecLike,
	Confirm: BitVecLike + std::ops::IndexMut<usize>,
	Might:   Mightsee<Filter, Confirm>
>(
	graph: &'a LeafGraph,
	work: &Work<Filter>,
	limiter: T,

	confirmsee: &mut Confirm,
	leaf: usize,
	mightsee: &Might
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
			mightsee.filter_for_unexplored(localsee, confirmsee)
		};

		if confirmsee[p] == true && any_unexplored == false {
			continue }

		if let Some(new_limiter) = limiter.traverse(graph, p) {
			//println!("   {{");
			confirmsee[p] = true;
			visited += recursive_leaf_flow(
				graph,
				&work,
				new_limiter,
				confirmsee,
				graph.portals[p].leaf_into,
				&*new_mightsee
			);
			//println!("   }}");
		}
	}
	return visited+1
}

pub fn flood(
	graph: &LeafGraph,
	portalfront: &Vec<BitVec>,
	portalflood: &mut BitVec,
	p: usize,
	leaf: usize)
{
	for rq in &graph.leaf_from[leaf] { // in portals out of leaf
		let q = *rq;
		if !portalfront[p][q] { continue }
		if portalflood[q] { continue }
		portalflood[q] = true;
		flood(graph, portalfront, portalflood, p, graph.portals[q].leaf_into);
	}
}


pub fn recursive_leaf_flow_simple<'a, T: Limiter<'a>>(
	graph: &'a LeafGraph,
	work: &Work<BitVec>,
	limiter: T,
	base: usize
) -> (BitVec, usize, usize, usize)
{
	let mut basevis: BitVec = Vec::new();
	basevis.resize_with(graph.portals.len(), Default::default);
	let mightsee: BitVec = {
		let guard = &epoch::pin();
		work.get_row(base, guard).clone()
	};
    let npapprox = bitvec_pop(&mightsee);
	assert_eq!(mightsee.len(), basevis.len());
	let nvisited = recursive_leaf_flow(
		&graph,
		&work,
		limiter,
		&mut basevis,
		graph.portals[base].leaf_into,
		&MaybeSparse::Dense(mightsee)
    );
    let nvisible = bitvec_pop(&basevis);
    (basevis, nvisible, nvisited-1, npapprox)
}

use std::time::Instant;


use crossbeam::thread;
use crossbeam::epoch::{self as epoch, Atomic, Owned, Guard};

pub struct Work<RowType> {
	approx: Vec<Atomic<RowType>>,
	first_visible:  AtomicU64,
	first_visited:  AtomicU64,
}

impl<R> Work<R> {
	fn get_row<'g>(&self, a: usize, guard: &'g Guard) -> &'g R {
		unsafe { self.approx[a].load(Ordering::SeqCst, guard).as_ref() }.expect("why is any row null?")
	}
	fn set_row(&self, a: usize, data: R) {
		self.approx[a].store(Owned::new(data), Ordering::SeqCst)
	}
}

fn winding_limiter_for<'a>(graph: &LeafGraph, p: usize, ond: Option<&'a Vec<(N, N)>>) -> WindingLimiter<'a> {
	let pt = &graph.portals[p];
	WindingLimiter(
		if let Some(nd) = ond {
			pt.plane.with_dist(nd[p].0)
		} else {
			pt.plane
		},
		pt.winding.clone(),
	 	None,
	 	ond)
}

fn distances_for(graph: &LeafGraph, i: usize) -> (N, N) {
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

	(
		if ahead_dist == std::f64::INFINITY { plane.3 } else { ahead_dist },
		if behind_dist == -std::f64::INFINITY { plane.3 } else { behind_dist }
	)
}

fn front_check(p: &Portal, q: &Portal) -> bool {
	!all_below(&q.winding.points, p.plane) &&
	!all_above(&p.winding.points, q.plane)
}

fn front_check_dist(p: &Portal, q: &Portal, d: N, e: N) -> bool {
	!all_below(&q.winding.points, p.plane.with_dist(d)) &&
	!all_above(&p.winding.points, q.plane.with_dist(e))
}

fn front_check_matrix(
	graph:        &LeafGraph,
	newdistances: &Vec<(N, N)>,
	portalfront:  &mut Vec<BitVec>,
	mode2:        bool
) {
	let np = graph.portals.len();
	for i in 0..np {
		for j in 0..np {
			portalfront[i][j] = i != j && if mode2 {
				front_check_dist(&graph.portals[i], &graph.portals[j],
					newdistances[i].0,
					newdistances[i].1)
			} else {
				front_check(&graph.portals[i], &graph.portals[j])
			}
		}
	}
}

use std::collections::BinaryHeap;

pub fn heap_merge(items: &Vec<&BitVec>) -> Vec<Vec<usize>> {
	let mut cchoices = Vec::new();
	let m = items[0].len();
	let mut h: BinaryHeap<(usize, usize, usize)> = BinaryHeap::new();
	for i in 0..items.len() {
		let ip = bitvec_pop(items[i]);
		if ip > 0 {
			h.push((m-ip, i, bitvec_pop(items[i])))
		}
	}

	while h.len() > 0 {
		let mut choices: Vec<usize> = Vec::new();
		let mut cbv = BitVec::new();
		cbv.resize(m, false);
		let mut cpop = 0;
		let mut c2pop = 0;

		while let Some((oldneguniques, index, itempop)) = h.pop() {
			let newuniques = bitvec_and_not_pop(items[index], &cbv);
			if oldneguniques != m-newuniques {
				h.push((m-newuniques, index, itempop));
				continue
			}
			let newpop = cpop + newuniques;
			let new2pop = c2pop + itempop;
			if choices.len() > 0 {
				let totalold = (cpop * cpop) + (c2pop * cpop) + 100 * itempop * itempop;
				let totalnew = (newpop * newpop) + (new2pop * newpop);
				if totalnew > totalold {
					break
				}
				println!("group {}, for {}th, add {:?}, score {} -> {} (with pop {}: {} -> {}, evals {} -> {})",
					cchoices.len(), choices.len(), index, totalold, totalnew, m, cpop, newpop, c2pop, new2pop);
			} else {
				println!("group {}, for {}th, add {:?}, initial",
					cchoices.len(), choices.len(), index);
			}
			choices.push(index);
			cbv = bitvec_or(&cbv, items[index]);
			cpop = newpop;
			c2pop = new2pop;
		}
		cchoices.push(choices);

		// reset all scores
		let mut h2 = BinaryHeap::new();
		for (_, i, ip) in h {
			h2.push((m-ip, i, ip));
		}
		h = h2;
	}

	cchoices
}

pub fn process_graph<'a>(graph: &'a LeafGraph) -> (Vec<BitVec>, Vec<BitVec>) {
	let nportals = graph.portals.len();

	println!("process_graph: init  nportals={}", nportals);

	let mut newdistances: Vec<(N, N)> = Vec::new();
	for i in 0..nportals {
		let dd = distances_for(graph, i);
		//println!("pplane +{:?}:-{:?}", dd.0-plane.3, plane.3-dd.1);
		newdistances.push(dd);
	}

	let mut portalflood: Vec<BitVec> = Vec::new();
	portalflood.resize_with(nportals, Default::default);

	let mut portalfront: Vec<BitVec> = Vec::new();
	portalfront.resize_with(nportals, Default::default);

	let mut portalvis: Vec<BitVec> = Vec::new();
	portalvis.resize_with(nportals, Default::default);

	for i in 0..nportals {
		portalfront[i].resize(nportals, false);
		portalflood[i].resize(nportals, false);
		portalvis[i].resize(nportals, false);
	}

	println!("process_graph: portalfront");
	let now = Instant::now();
	front_check_matrix(&graph, &newdistances, &mut portalfront, false);
	println!("process_graph: portalfront took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);


	println!("process_graph: portalflood");
	let now = Instant::now();
	for i in 0..nportals {
		flood(&graph, &portalfront, &mut portalflood[i], i, graph.portals[i].leaf_into);
	}
	println!("process_graph: portalflood took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);

	let workunit = move |graph: &LeafGraph, work: &Work<BitVec>, p: usize| {
		let limiter_w = winding_limiter_for(&*graph, p, None); // Some(&newdistances)
		// let limiter_fm = FMLimiter(vec![&graph.portals[p].winding]);
		// let limiter =  (limiter_w, limiter_fm);

		// let now = Instant::now();
    	let (row, visible_w, visited_w, cc) = recursive_leaf_flow_simple(
    		&graph, &work, limiter_w.clone(), p);
    	// let (_, visited_w2) = recursive_leaf_flow_simple(
    	// 	&graph, &woget_rowrk, limiter_w2, p);
    	// let (row, visited) = recursive_leaf_flow_simple(
    	// 	&graph, &work, limiter, p);
    	// let (_, visited2) = recursive_leaf_flow_simple(
    	// 	&graph, &work, limiter2, p);

    	work.set_row(p, row);

    	work.first_visible.fetch_add(visible_w as u64, Ordering::SeqCst);
    	work.first_visited.fetch_add(visited_w as u64, Ordering::SeqCst);

		// if visited_w != visible_w {
		// 	println!(" (visited {:?}/{:?} leafs, of which {:?} are visible)", visited_w, cc, visible_w);
		// } else {
		// 	println!(" (visited {:?}/{:?} leafs)", visited_w, cc);
		// }
    };
	// for p in 0..nportals {
	// 	workunit(p);
	// }

	//let groups: Vec<Vec<usize>> = heap_merge(&portalflood.iter().map(|a|a).collect());
	let mut nonzeroes: Vec<usize> = Vec::new();
	for i in 0..nportals {
		if graph.leaf_from[graph.portals[i].leaf_into].len() > 0 {
			nonzeroes.push(i);
		}
	}
	let groups: Vec<Vec<usize>> = vec![nonzeroes];

	for (i, g) in groups.iter().enumerate() {
		let mut sub_indices_set = HashSet::<usize>::new();
		for p in g {
			for q in 0..portalflood[*p].len() {
				if portalflood[*p][q] {
					sub_indices_set.insert(q);
				}
			}
		}
		let sub_indices: Vec<usize> = sub_indices_set.into_iter().collect();
		let nlportals = sub_indices.len();
		println!("subgraph {}: entry through {}, total {}", i, g.len(), nlportals);
		let sub_graph = renumber(graph, &sub_indices);
		let sub_portalflood = reduce_matrix(&portalflood, &sub_indices);
		let sub_portalvis: Vec<Atomic<BitVec>> =
			sub_portalflood.iter().map(|row| Atomic::new(row.clone())).collect();

		let work = Work {
			approx: sub_portalvis,
			first_visible: 0.into(),
			first_visited: 0.into(),
		};

		let now = Instant::now();
		thread::scope(|s| {
			let rwork = &work;
			let rworkunit = &workunit;
			let rsubgraph = &sub_graph;
			for i in 0..4 {
				s.spawn(move |_| {
				 	for p in 0..nlportals {
				 		if p % 4 == i {
							// print!("recursive_leaf_flow {:?}/{:?}\n", p, nlportals);
							rworkunit(rsubgraph, rwork, p);
						}
					}
				});
			}
		}).unwrap();
		println!("process_graph: subgraph took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);

		let owned_rows: Vec<Box<BitVec>> = work.approx.into_iter().map(|atomic|
			unsafe { atomic.load(Ordering::Acquire, &epoch::pin()).into_owned().into_box() }).collect();

		let owned_rows = owned_rows.into_iter().map(|b|*b).collect();

		reinject_matrix(&mut portalvis, owned_rows, &sub_indices);

		println!("{} {}",
			work.first_visible.load(Ordering::Acquire),
			work.first_visited.load(Ordering::Acquire)
		);
	}

	(portalvis_to_leafvis(&graph, &portalflood), portalvis_to_leafvis(&graph, &portalvis))
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
