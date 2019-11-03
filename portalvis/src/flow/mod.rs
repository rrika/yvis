pub mod limiter;
pub mod mightsee;

use crate::bits::BitVecLike;
use crate::bits::chunked_convert;
use crate::bits::ChunkedBits;
use crate::bits::DenseBits;
use crate::bits::SetBit;
use crate::flow::limiter::*;
use crate::flow::mightsee::Mightsee;
use crate::flow::mightsee::MightseeResult;
use crate::geometry::all_above;
use crate::geometry::all_below;
use crate::geometry::N;
use crate::geometry::ndot;
use crate::graph::LeafGraph;
use crate::graph::Portal;
use crate::report;
use std::collections::HashSet;
use std::sync::atomic::AtomicU64;
use std::sync::atomic::Ordering;
use std::time::Instant;
use crossbeam::thread;
use crossbeam::epoch::{self as epoch, Atomic, Owned, Guard};

// ================================================================ //

pub struct Work<RowType> {
	approx: Vec<Atomic<RowType>>,
	// visible:  AtomicU64,
	// wasted:  AtomicU64,
	// visited:  AtomicU64,
}

impl<R> Work<R> {
	fn get_row<'g>(&self, a: usize, guard: &'g Guard) -> &'g R {
		unsafe { self.approx[a].load(Ordering::SeqCst, guard).as_ref() }.expect("why is any row null?")
	}
	fn set_row(&self, a: usize, data: R) {
		self.approx[a].store(Owned::new(data), Ordering::SeqCst)
	}
	fn report_traversal(&self, t: &report::Traversal) {
		// println!("{:?}", t);
		// self.visible.fetch_add(t.result as u64, Ordering::SeqCst);
		// self.wasted.fetch_add(t.wasted as u64, Ordering::SeqCst);
		// self.visited.fetch_add(t.steps as u64, Ordering::SeqCst);
	}
}

pub fn reduce_matrix(m: &Vec<DenseBits>, r: &Vec<usize>) -> Vec<DenseBits> {
	r.iter().map(|x| {
		let mut row = DenseBits(Vec::new());
		row.0.resize((r.len()+63)/64, 0u64);
		for (i, y) in r.iter().enumerate() {
			if m[*x][*y] {
				row.set(i)
			}
		}
		row
	}).collect()
}

pub fn reinject_matrix(m: &mut Vec<DenseBits>, n: Vec<DenseBits>, p: &Vec<usize>) {
	for i in 0..p.len() {
		for j in 0..p.len() {
			if n[i][j] {
				m[p[i]].set(p[j])
			}
		}
	}
}

pub fn matrix_make_symmetric(matrix: &mut Vec<DenseBits>) {
	let nportals = matrix.len();
	for p in 1..nportals {
		for q in 0..p-1 {
			let a = matrix[q^1][p^1];
			let b = matrix[p][q];
			if a != b {
				matrix[q^1].clear(p^1);
				matrix[p].clear(q);
			}
		}
	}
}

// ================================================================ //

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
	//let r = !(p.plane == -q.plane);
	let s = !all_below(&q.winding.points, p.plane);
	let t = !all_above(&p.winding.points, q.plane);
	//println!("= {} = {} && {} && {} ", r&&s&&t, r, s, t);
	//r&&s&&t
	s&&t
}

fn front_check_dist(p: &Portal, q: &Portal, d: N, e: N) -> bool {
	!all_below(&q.winding.points, p.plane.with_dist(d)) &&
	!all_above(&p.winding.points, q.plane.with_dist(e))
}

pub fn calc_portalfront(
	graph:        &LeafGraph,
	portalfront:  &mut Vec<DenseBits>,
	newdistances: Option<&Vec<(N, N)>>
) {
	let nportals = graph.portals.len();
	for i in 0..nportals {
		for j in 0..nportals {
			//print!("front check {} {} ", i, j);
			let f = i != j && if let Some(newdistances) = newdistances {
				front_check_dist(&graph.portals[i], &graph.portals[j],
					newdistances[i].0,
					newdistances[i].1)
			} else {
				front_check(&graph.portals[i], &graph.portals[j])
			};
			if f {
				portalfront[i].set(j)
			}
		}
	}
}

// ================================================================ //

pub fn flood(
	graph: &LeafGraph,
	portalfront: &Vec<DenseBits>,
	portalflood: &mut DenseBits,
	p: usize,
	leaf: usize)
{
	for rq in &graph.leaf_from[leaf] { // in portals out of leaf
		let q = *rq;
		let nleaf = graph.portals[q].leaf_into;
		//assert_eq!(portalfront[q^1][p^1], portalfront[p][q]);
		if !portalfront[p][q] { continue }
		if portalflood[q] { continue }
		portalflood.set(q);
		flood(graph, portalfront, portalflood, p, nleaf);
	}
}

pub fn calc_portalflood(
	graph:        &LeafGraph,
	portalflood:  &mut Vec<DenseBits>,
	portalfront:  &Vec<DenseBits>
) {
	let nportals = graph.portals.len();
	for i in 0..nportals {
		flood(&graph, &portalfront, &mut portalflood[i], i, graph.portals[i].leaf_into);
	}
}

// ================================================================ //

pub fn recursive_leaf_flow<
	'a,
	T: Limiter<'a> + Clone,
	Filter:  BitVecLike,
	Confirm: BitVecLike + SetBit,
	//Might:   Mightsee<Filter, Confirm>
>(
	graph: &'a LeafGraph,
	work: &Work<Filter>,
	limiter: T,

	confirmsee: &mut Confirm,
	leaf: usize,
	//mightsee: &Might,
	mightsee: &ChunkedBits,
	depth: usize
) -> (usize, usize, bool) where ChunkedBits: Mightsee<Filter, Confirm> {
	let mut any_confirmed = false;
	let mut visited = 0;
	let mut wasted = 0;
	for pp in &graph.leaf_from[leaf] { // in portals out of leaf
		let p = *pp;
		//println!("  {:?}", p);

		if mightsee[p] == false {
			//println!("    may not see");
			continue }

		//let local_mightsee: Might;
		let local_mightsee: ChunkedBits;

		let now = Instant::now();
		let rnms: Option<&ChunkedBits> = match {
			let guard = &epoch::pin();
			let localsee = work.get_row(p, guard);
			mightsee.filter_for_unexplored(localsee, confirmsee)
		} {
			MightseeResult::None => {
				if confirmsee[p] == true { continue }
				None
			},
			MightseeResult::Unchanged => {
				Some(&mightsee)
			},
			MightseeResult::Reduced(r) => {
				local_mightsee = r;
				Some(&local_mightsee)
			}
		};
		let mightsee_time = now.elapsed();

		if let Some(new_limiter) = limiter.traverse(graph, p) {
			let limiter_time = now.elapsed() - mightsee_time;
			//println!("   {{");
			if !confirmsee[p] {
				any_confirmed = true;
			}
			confirmsee.set(p);
			if let Some(new_mightsee) = rnms {
				let (sub_visited, sub_wasted, sub_anyconfirmed) = recursive_leaf_flow(
					graph,
					&work,
					new_limiter,
					confirmsee,
					graph.portals[p].leaf_into,
					new_mightsee,
					depth + 1
				);
				visited += sub_visited;
				wasted += sub_wasted;
				any_confirmed |= sub_anyconfirmed;
			}
			//println!("   }}");

			// print!("depth={} mightsee={} {}ns {}ns",
			// 	depth,
			// 	mightsee.0.len()*64,
			// 	mightsee_time.as_nanos() as u64,
			// 	limiter_time.as_nanos() as u64);
			// limiter.stats();
		}
	}
	return (
		visited+1,
		if any_confirmed { wasted } else { wasted+1 },
		any_confirmed
	)
}

pub fn recursive_leaf_flow_simple<'a, T: Limiter<'a>+Clone, F: Fn(usize)->T>(
	graph: &'a LeafGraph,
	work: &Work<DenseBits>,
	//limiter: T,
	glimiter: F,
	bases: &[usize]
) -> (DenseBits, report::Traversal)
{
	let mut npapprox = 0;
	let mut nvisited = 0;
	let mut nwasted = 0;
	let nportals = graph.portals.len();
	let mut basevis: DenseBits = DenseBits(Vec::new());
	basevis.0.resize((nportals + 63)/64, 0); // nn
	for base in bases {
		let mightsee: DenseBits = {
			let guard = &epoch::pin();
			work.get_row(*base, guard).clone()
		};
		npapprox = mightsee.count_ones();
		let (snvisited, snwasted, _sanyconf) = recursive_leaf_flow(
			&graph,
			&work,
			glimiter(*base),
			&mut basevis,
			graph.portals[*base].leaf_into,
			//&sparse_convert(&mightsee, nportals),
			&chunked_convert(&mightsee),
			//&mightsee,
			0
		);
		nvisited += snvisited;
		nwasted += snwasted;
	}
	let report = report::Traversal {
		approx:  npapprox,
		steps:   nvisited-1,
		result:  basevis.count_ones(),
		wasted:  nwasted,
		depth:   0
	};
	(basevis, report)
}

fn winding_limiter_for<'a>(
	graph: &LeafGraph,
	p:     usize,
	ond:   Option<&'a Vec<(N, N)>>) -> WindingLimiter<'a>
{
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

fn workunit(graph: &LeafGraph, work: &Work<DenseBits>, p: usize) {
	let pt = &graph.portals[p];
	let _limiter_w = winding_limiter_for(graph, p, None); // Some(&newdistances)
	let _limiter_chop = ChopLimiter::<crate::axisgeom::Cuboid>(pt.plane, (&pt.winding).into(), None);
	// let limiter_fm = FMLimiter(vec![&graph.portals[p].winding]);
	// let limiter =  (limiter_w, limiter_fm);
	// let limiter = CheckLimiter(_limiter_chop, Some(_limiter_w));
	let limiter = _limiter_chop;
	// let limiter = _limiter_w;

	let (row, traversal) = recursive_leaf_flow_simple(
		&graph, &work, |_| limiter.clone(), &[p]);
	// let (_, visited_w2) = recursive_leaf_flow_simple(
	// 	&graph, &work, limiter_w2, p);
	// let (row, visited) = recursive_leaf_flow_simple(
	// 	&graph, &work, limiter, p);
	// let (_, visited2) = recursive_leaf_flow_simple(
	// 	&graph, &work, limiter2, p);
	work.set_row(p, row);
	work.report_traversal(&traversal);

}

pub fn calc_portalvis(
	graph:        &LeafGraph,
	portalvis:    &mut Vec<DenseBits>,
	portalflood:  &Vec<DenseBits>
) {
	let nportals = graph.portals.len();

	// for p in 0..nportals {
	// 	workunit(p);
	// }

	// let groups: Vec<Vec<usize>> = heap_merge(&portalflood.iter().map(|a|a).collect(), nportals);
	let mut nonzeroes: Vec<usize> = Vec::new();
	for i in 0..nportals {
		// if graph.leaf_from[graph.portals[i].leaf_into].len() > 0 {
			nonzeroes.push(i);
		// }
	}
	let groups: Vec<Vec<usize>> = vec![nonzeroes];

	for i in 0..nportals {
		portalvis[i].0 = portalflood[i].0.clone();
	}

	for (i, g) in groups.iter().enumerate() {
		let mut sub_indices_set: HashSet<usize> = g.clone().into_iter().collect::<HashSet<usize>>();
		for p in g {
			for q in 0..nportals {
				if portalflood[*p][q] {
					sub_indices_set.insert(q);
				}
			}
		}
		for p in g {
			sub_indices_set.remove(p);
		}
		assert_eq!(sub_indices_set.len(), 0);
		let mut sub_indices: Vec<usize> = g.clone();
		let nlportals = sub_indices.len();
		sub_indices.extend(sub_indices_set.into_iter().collect::<Vec<usize>>());

		println!("process_graph: subgraph {}: entry through {}, total {}", i, g.len(), nlportals);
		// let sub_graph = renumber(graph, &sub_indices);
		// let sub_local_indices = (0..nlportals).collect::<Vec<usize>>();
		// let sub_portalvis_ = reduce_matrix(&portalvis, &sub_indices);
		let sub_graph = &graph;
		let sub_local_indices = &sub_indices;
		let sub_portalvis_ = &portalvis;

		let sub_portalvis: Vec<Atomic<DenseBits>> =
			sub_portalvis_.iter().map(|row| Atomic::new(row.clone())).collect();

		let work = Work {
			approx: sub_portalvis,
			// visible: 0.into(),
			// wasted: 0.into(),
			// visited: 0.into(),
		};

		// println!("subgraph {}: pre pass", i);
		// let (row, visible_w, visited_w, cc) = recursive_leaf_flow_simple(
		// 	&graph, &work, |p|winding_limiter_for(&*graph, p, None),
		// 	sub_local_indices);

		let nthreads = 4;

		println!("process_graph: subgraph {}: main pass", i);
		let now = Instant::now();
		thread::scope(|s| {
			let rwork = &work;
			let rsubgraph = &sub_graph;
			for i in 0..nthreads {
				s.spawn(move |_| {
					for p in 0..sub_local_indices.len() {
						if p % nthreads == i {
							// print!("recursive_leaf_flow {:?}/{:?}\n", p, nlportals);
							workunit(rsubgraph, rwork, sub_local_indices[p]);
						}
					}
				});
			}
		}).unwrap();
		println!("process_graph: subgraph took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);

		let owned_rows: Vec<Box<DenseBits>> = work.approx.into_iter().map(|atomic|
			unsafe { atomic.load(Ordering::Acquire, &epoch::pin()).into_owned().into_box() }).collect();

		let mut owned_rows: Vec<DenseBits> = owned_rows.into_iter().map(|b|*b).collect();

		// reinject_matrix(&mut portalvis, owned_rows, &sub_indices);
		for p in sub_indices {
			std::mem::swap(&mut portalvis[p], &mut owned_rows[p]);
		}

		// println!("process_graph: stats: visible={} wasted={} visited={}",
		// 	work.visible.load(Ordering::Acquire),
		// 	work.wasted.load(Ordering::Acquire),
		// 	work.visited.load(Ordering::Acquire)
		// );
	}
}

// ================================================================ //

pub fn process_graph<'a>(graph: &'a LeafGraph) -> (Vec<DenseBits>, Vec<DenseBits>) {
	let nportals = graph.portals.len();

	println!("process_graph: init  nportals={}", nportals);

	// let mut newdistances: Vec<(N, N)> = Vec::new();
	// for i in 0..nportals {
	// 	let dd = distances_for(graph, i);
	// 	newdistances.push(dd);
	// }

	let mkbitvec = |n| { let mut v = Vec::new(); v.resize((n+63)/64, 0u64); DenseBits(v) };

	let mut portalflood: Vec<DenseBits> = graph.square(mkbitvec);
	let mut portalfront: Vec<DenseBits> = graph.square(mkbitvec);
	let mut portalvis: Vec<DenseBits> = graph.square(mkbitvec);
	let mut portalvis2: Vec<DenseBits> = graph.square(mkbitvec);

	println!("process_graph: portalfront");
	let now = Instant::now();
	calc_portalfront(&graph, &mut portalfront, None); //Some(&newdistances));
	println!("process_graph: portalfront took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);


	println!("process_graph: portalflood");
	let now = Instant::now();
	calc_portalflood(&graph, &mut portalflood, &portalfront);
	matrix_make_symmetric(&mut portalflood);
	println!("process_graph: portalflood took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);

	println!("process_graph: portalvis");
	let now = Instant::now();
	calc_portalvis(&graph, &mut portalvis, &portalflood);
	matrix_make_symmetric(&mut portalvis);
	calc_portalvis(&graph, &mut portalvis2, &portalvis);
	matrix_make_symmetric(&mut portalvis2);
	// matrix_make_symmetric(&mut portalvis);
	println!("process_graph: portalvis took {}ms", now.elapsed().as_micros() as u64 as f64 / 1000.0);

	(portalflood, portalvis2)
}
