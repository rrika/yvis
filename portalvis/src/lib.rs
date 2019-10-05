use std::fmt::Debug;
use std::sync::atomic::AtomicBool;
use std::sync::atomic::Ordering;

type N = f64;
type Row = [N; 6];

#[derive(Copy, Clone, Debug)]
pub enum Relation {
	GeZero,
	GtZero
}

type System<T> = Vec<(Row, Relation, T)>;

pub struct ReconstructionInfo<T: DeriviationTracker> {
	gt: Vec<(Row, T::T)>,
	ge: Vec<(Row, T::T)>,
	lt: Vec<(Row, T::T)>,
	le: Vec<(Row, T::T)>,
}

pub fn dot(lhs: Row, rhs: Row) -> N {
	lhs[0] * rhs[0] + 
	lhs[1] * rhs[1] + 
	lhs[2] * rhs[2] + 
	lhs[3] * rhs[3] + 
	lhs[4] * rhs[4] + 
	lhs[5] * rhs[5]
}

pub fn reconstruct<T: DeriviationTracker>(lhs: Row, recon: &ReconstructionInfo<T>, debugname: &str, deriv: &mut T)
	-> Result<N, T::T>
{
	let gt: Option<(N, T::T)> = recon.gt.iter().map(|(rhs, t)| (dot(lhs, *rhs), *t)).min_by(|a, b| (a.0).partial_cmp(&b.0).unwrap());
	let ge: Option<(N, T::T)> = recon.ge.iter().map(|(rhs, t)| (dot(lhs, *rhs), *t)).min_by(|a, b| (a.0).partial_cmp(&b.0).unwrap());
	let lt: Option<(N, T::T)> = recon.lt.iter().map(|(rhs, t)| (dot(lhs, *rhs), *t)).max_by(|a, b| (a.0).partial_cmp(&b.0).unwrap());
	let le: Option<(N, T::T)> = recon.le.iter().map(|(rhs, t)| (dot(lhs, *rhs), *t)).max_by(|a, b| (a.0).partial_cmp(&b.0).unwrap());

	let g = match (gt, ge) {
		(None,    None)                  => None,
		(Some(v), None)                  => Some((false, v)),
		(None,    Some(w))               => Some((true, w)),
		(Some(v), Some(w)) if v.0 <= w.0 => Some((false, v)),
		(Some(v), Some(w)) if v.0 >  w.0 => Some((true, w)),
		_ => panic!("NaN somewhere")
	};
	let l = match (lt, le) {
		(None,    None)                  => None,
		(Some(v), None)                  => Some((false, v)),
		(None,    Some(w))               => Some((true, w)),
		(Some(v), Some(w)) if v.0 >= w.0 => Some((false, v)),
		(Some(v), Some(w)) if v.0 <  w.0 => Some((true, w)),
		_ => panic!("NaN somewhere")
	};
	println!("{:?} ... {} ... {:?}", l, debugname, g);
	match (g, l) {
		(None, None)                                     => Ok(0.0),
		(None, Some((eq, w)))                            => Ok(if eq { w.0 } else { w.0 + 1.0 }),
		(Some((eq, v)), None)                            => Ok(if eq { v.0 } else { v.0 - 1.0 }),
		(Some((true, v)), Some((true, w))) if v.0 >= w.0 => Ok(0.5 * (v.0+w.0)),
		(Some((_, v)), Some((_, w)))       if v.0 >  w.0 => Ok(0.5 * (v.0+w.0)),
		(Some((_, v)), Some((_, w)))                     => Err(deriv.combine(v.1, w.1)),
	}
}

pub fn divide_row(lhs: Row, rhs: N) -> Row {
	[
		lhs[0]/rhs,
		lhs[1]/rhs,
		lhs[2]/rhs,
		lhs[3]/rhs,
		lhs[4]/rhs,
		lhs[5]/rhs
	]
}

pub fn subtract_row(lhs: Row, rhs: Row) -> Row {
	[
		lhs[0] - rhs[0],
		lhs[1] - rhs[1],
		lhs[2] - rhs[2],
		lhs[3] - rhs[3],
		lhs[4] - rhs[4],
		lhs[5] - rhs[5]
	]
}

pub trait DeriviationTracker {
	type T: Copy + Clone + Debug;
	fn new_item(&mut self) -> Self::T;
	fn combine(&mut self, a: Self::T, b: Self::T) -> Self::T;
}

pub struct BitfieldTracker(usize);
impl DeriviationTracker for BitfieldTracker {
	type T = u64;
	fn new_item(&mut self) -> u64 {
		self.0 += 1;
		1u64 << self.0-1
	}
	fn combine(&mut self, a: u64, b: u64) -> u64 { a | b }
}

pub fn any_solution<T: DeriviationTracker>(
	ineqs: &System<T::T>,
	first_axis: usize,
	last_axis: usize,
	deriv: &mut T) -> Result<Row, T::T>
{
	let mut ineqs = ineqs.clone();
	let mut reconstruction: Vec<(usize, ReconstructionInfo<T>)> = Vec::new();

	for axis in first_axis..=last_axis {
		for ineq in &ineqs { println!("{:?}", ineq); } println!("--");
		reconstruction.push((axis, {
			let (recon, new_ineqs) = fourier_motzkin_elimination(&ineqs, axis, deriv);
			ineqs = new_ineqs;
			recon
		}));
	}
	let mut solution: Row = [0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
	for (axis, rinfo) in reconstruction.iter().rev() {
		solution[*axis] = reconstruct::<T>(solution, &rinfo, &["x", "y"][*axis], deriv)?;
		println!("sol {:?}", solution);
	}
	Ok(solution)
}

type BitVec = Vec<bool>;


#[derive(Clone, Debug)]
pub struct Winding { points: Vec<[N; 3]> }

#[derive(Copy, Clone, Debug, Default)]
pub struct Plane(N, N, N, N);

impl std::cmp::PartialEq for Plane {
	fn eq(&self, rhs: &Plane) -> bool {
		self.0 * rhs.3 == self.3 * rhs.0 && 
		self.1 * rhs.3 == self.3 * rhs.1 && 
		self.2 * rhs.3 == self.3 * rhs.2
	}
}

impl std::ops::Neg for Plane {
	type Output = Plane;
	fn neg(self) -> Self::Output { Plane(-self.0, -self.1, -self.2, -self.3) }
}

impl std::ops::Mul<[N; 3]> for Plane {
	type Output = N;
	fn mul(self, rhs: [N; 3]) -> Self::Output {
		  self.0 * rhs[0] 
		+ self.1 * rhs[1]
		+ self.2 * rhs[2]
		- self.3
	}
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

pub fn chop_winding(w: &Winding, p: Plane) -> Option<Winding> {
	//println!("chop_winding {:?} {:?}", p, w);
	let mut nw: Vec<[N; 3]> = Vec::new();
	let mut add_point = |p: [N; 3]| {
		if nw.len() >= 2 {
			let n = nw[nw.len()-2];
			let o = nw[nw.len()-1];
			if ncross3(n, o, p) == [0.0, 0.0, 0.0] {
				let last = nw.len()-1;
				nw[last] = p;
				return
			}
		}
		nw.push(p);
	};
	let mut i = w.points.len()-1;
	for j in 0..w.points.len() {
		let c = w.points[i];
		let d = w.points[j];
		i = j;
		let a = p * c;
		let b = p * d;
		if a < 0.0 && b < 0.0 {
			continue
		} else if a >= 0.0 && b >= 0.0 {
			add_point(d)
		} else {

			let f = a / (a-b);
			let mid = nadd(c, nmul(nsub(d, c), f));

			if a < 0.0 && b >= 0.0 {
				add_point(mid);
				add_point(d);
				
			} else if a >= 0.0 && b < 0.0 {
				add_point(mid);
			} else {
				unimplemented!();
			}
		}
	}
	//println!("  chop -> {:?}", nw);
	if nw.len() >= 3 {
		Some(Winding{points: nw})
	} else {
		None
	}
}

pub fn ncross(a: [N; 3], b: [N; 3]) -> [N; 3] { [
	a[1]*b[2]-b[1]*a[2],
	a[2]*b[0]-b[2]*a[0],
	a[0]*b[1]-b[0]*a[1]
] }
pub fn ncross3(a: [N; 3], b: [N; 3], c: [N; 3]) -> [N; 3] { ncross(nsub(a, b), nsub(c, b)) }
pub fn nadd(a: [N; 3], b: [N; 3]) -> [N; 3] { [a[0]+b[0], a[1]+b[1], a[2]+b[2]] }
pub fn nsub(a: [N; 3], b: [N; 3]) -> [N; 3] { [a[0]-b[0], a[1]-b[1], a[2]-b[2]] }
pub fn nmul(a: [N; 3], b: N) -> [N; 3] { [a[0]*b, a[1]*b, a[2]*b]	}
pub fn ndot(a: [N; 3], b: [N; 3]) -> N { a[0]*b[0] + a[1]*b[1] + a[2]*b[2] }
pub fn plane_through(a: [N; 3], b: [N; 3], c: [N; 3]) -> Plane {
	let n = ncross3(a, b, c);
	let d = ndot(n, a);
	Plane(n[0], n[1], n[2], d)
}

pub fn all_above(points: &[[N; 3]], plane: Plane) -> bool {
	for p in points { if plane * *p < 0.0 { return false } }
	true
}
pub fn all_below(points: &[[N; 3]], plane: Plane) -> bool {
	for p in points { if plane * *p > 0.0 { return false } }
	true
}

pub fn clip_hourglass(source: &Winding, pass: &Winding, target: &Winding) -> Option<Winding> {
	let mut target_clipped = target.clone();

	for i in 0..source.points.len() {
		let j = (i+1) % source.points.len();
		for k in 0..pass.points.len() {
			let l = (k+1) % pass.points.len();

			let p1 = plane_through(
				source.points[i],
				source.points[j],
				pass.points[k]);

			let p2 = plane_through(
				source.points[i],
				pass.points[k],
				pass.points[l]);


			// all_below(&source.points, p1) is true
			if all_above(&pass.points, p1) {
				target_clipped = chop_winding(&target_clipped, p1)?;
			}

			// all_above(&pass.points, p2) is true
			if all_below(&source.points, p2) {
				target_clipped = chop_winding(&target_clipped, p2)?;
			}
		}
	}
	Some(target_clipped)
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

trait Limiter {
	fn traverse(&self, graph: &LeafGraph, portal: usize) -> Option<Self> where Self: Sized;
}

impl Limiter for MightSeeLimiter {
	fn traverse(&self, _graph: &LeafGraph, _portal: usize) -> Option<MightSeeLimiter> {
		panic!()
	}
}

impl Limiter for WindingLimiter {
	fn traverse(&self, graph: &LeafGraph, portal: usize) -> Option<Self> {
		let (source_plane, mut source) = (self.0, self.1.clone());
		let (target_plane, mut target) = (graph.plane[portal], graph.winding[portal].clone());
		//println!("   target_plane = {:?}", target_plane);
		if source_plane == -target_plane {
			//println!("   u-turn");
			return None
		}
		target = chop_winding(&target, -source_plane)?;
		source = chop_winding(&source, target_plane)?;
		if let Some(pass) = &self.2 {
			target = clip_hourglass(&source, &pass, &target)?;
		}
		Some(WindingLimiter(source_plane, source, Some(target)))
	}
}

pub struct Watcher(u8, u8);
pub struct LearningLimiter(Vec<Watcher>, BitVec);


impl Limiter for LearningLimiter {
	fn traverse(&self, _graph: &LeafGraph, _portal: usize) -> Option<Self> {
		panic!()
	}
}

pub fn recursive_leaf_flow(
	graph: &LeafGraph,
	base: usize,
	leaf: usize,
	mightsee: &BitVec,
	winding_limiter: WindingLimiter,

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
	(visited + 1)
}

pub fn fourier_motzkin_elimination<T: DeriviationTracker>(
	ineqs: &System<T::T>,
	axis: usize,
	deriv: &mut T
) -> (ReconstructionInfo<T>, System<T::T>) {
	let mut new_ineqs: System<T::T> = Vec::new();
	let mut reconstruct = ReconstructionInfo {
		gt: Vec::new(),
		ge: Vec::new(),
		lt: Vec::new(),
		le: Vec::new(),
	};

	for (row, relation, deriv) in ineqs {
		let v = row[axis];
		if v == 0.into() {
			new_ineqs.push((*row, *relation, *deriv));
			continue
		}
		let ar = divide_row(*row, -v);
		match (v > 0.into(), relation) {
			(true,  Relation::GtZero) => reconstruct.lt.push((ar, *deriv)),
			(false, Relation::GtZero) => reconstruct.gt.push((ar, *deriv)),
			(true,  Relation::GeZero) => reconstruct.le.push((ar, *deriv)),
			(false, Relation::GeZero) => reconstruct.ge.push((ar, *deriv)),
		}

	}
	for (above, below, rel) in [
		(&reconstruct.gt, &reconstruct.lt, Relation::GtZero),
		(&reconstruct.gt, &reconstruct.le, Relation::GtZero),
		(&reconstruct.ge, &reconstruct.lt, Relation::GtZero),
		(&reconstruct.ge, &reconstruct.le, Relation::GeZero)
	].iter() {
		for (oneabove, at) in *above {
			for (onebelow, bt) in *below {
				let mut new_row = subtract_row(*oneabove, *onebelow);
				new_row[axis] = 0.0;
				new_ineqs.push((new_row, *rel, deriv.combine(*at, *bt)));
			}
		}
	}
	(reconstruct, new_ineqs)
}

// let mut graph: LeafGraph = Default::default();

#[derive(Clone, Debug, Default)]
pub struct PRTLine(usize, usize, Vec<[f32; 3]>);

pub fn parse_prt(data: &str) -> Result<(usize, usize, Vec<PRTLine>), &'static str> {

	let mut prtlines = Vec::<PRTLine>::new();

	let mut line_iter = data.split("\r\n");
	let magic = line_iter.next().ok_or("prt line 1 empty")?;
	let nleafs = line_iter.next().ok_or("prt line 2 empty")?;
	let nportals = line_iter.next().ok_or("prt line 3 empty")?;

	if magic != "PRT1" { return Err("first line isn't PRT1") }
	let nleafs_i = nleafs.parse::<usize>().map_err(|_| "nleafs not an integer")?;
	let nportals_i = nportals.parse::<usize>().map_err(|_| "nportals not an integer")?;

	let mut nleafs_implicit = 0;

	for _p in 0..nportals_i {
		let line = line_iter.next().ok_or("file truncated")?;
		let mut parts = line.split("(");
		let mut nlab = parts.next().ok_or("prt line without points")?.split(" ");
		let npoints = nlab.next().ok_or("prt line missing npoints")?;
		let leaf_a = nlab.next().ok_or("prt line missing leaf_a")?;
		let leaf_b = nlab.next().ok_or("prt line missing leaf_b")?;

		let mut winding = Vec::<[f32; 3]>::new();

		for _i in 0..npoints.parse::<usize>().map_err(|_| "npoints not an integer")? {
			// expected format is (x y z ) with space after the z
			let mut xyz = parts.next().ok_or("prt line with fewer points than declared")?.split(" ");
			let x = xyz.next().ok_or("prt point lacking x")?;
			let y = xyz.next().ok_or("prt point lacking y")?;
			let z = xyz.next().ok_or("prt point lacking z")?;

			winding.push([
				x.parse().map_err(|_| "prt point x not an integer")?,
				y.parse().map_err(|_| "prt point y not an integer")?,
				z.parse().map_err(|_| "prt point z not an integer")?
			]);
		}

		let leaf_a_i = leaf_a.parse::<usize>().map_err(|_| "prt leaf_a not an integer")?;
		let leaf_b_i = leaf_b.parse::<usize>().map_err(|_| "prt leaf_b not an integer")?;

		prtlines.push(PRTLine(leaf_a_i, leaf_b_i, winding));

		nleafs_implicit = std::cmp::max(nleafs_implicit, leaf_a_i+1);
		nleafs_implicit = std::cmp::max(nleafs_implicit, leaf_b_i+1);
	}

	Ok((nleafs_i, nleafs_implicit, prtlines))
}

pub fn prtlines_to_graph(nleafs: usize, lines: &Vec<PRTLine>) -> LeafGraph {
	let mut graph = LeafGraph::default();

	let nportals = 2 * lines.len(); // prt portals are two-sided
	graph.leaf_from.reserve(nleafs);
	graph.leaf_into.reserve(nportals);
	graph.plane.reserve(nportals);
	graph.winding.reserve(nportals);
	for _i in 0..nleafs { graph.leaf_from.push(Vec::new()) }

	for PRTLine(src, dst, points) in lines {
		let mut points: Vec<[N; 3]> = points.iter().map(|a|[a[0] as N, a[1] as N, a[2] as N]).collect();

		graph.leaf_from[*dst].push(graph.leaf_into.len()); graph.leaf_into.push(*src);
		graph.leaf_from[*src].push(graph.leaf_into.len()); graph.leaf_into.push(*dst);
		let plane = plane_through(points[0], points[1], points[2]);
		graph.winding.push(Winding{points: points.clone()}); points.reverse();
		graph.winding.push(Winding{points: points});
		graph.plane.push(-plane);
		graph.plane.push(plane);
	}

	graph
}

pub fn process_graph(graph: &LeafGraph) {
	let nportals2 = graph.winding.len();

	let mut portalflood: Vec<BitVec> = Vec::new();
	let mut portalvis: Vec<BitVec> = Vec::new();
	let mut progress: Vec<std::sync::atomic::AtomicBool> = Vec::new();
	portalflood.resize_with(nportals2, Default::default);
	portalvis.resize_with(nportals2, Default::default);
	progress.resize_with(nportals2, || AtomicBool::new(false));

	for i in 0..nportals2 {
		portalflood[i].resize_with(nportals2, || true);
		portalvis[i].resize_with(nportals2, Default::default);
	}


	let mut visited = 0;
	for p in 0..nportals2 {
		println!("recursive_leaf_flow {:?}/{:?} (visited {:?} leafs)", p, nportals2, visited);
		let mightsee = &portalflood[p];
		let winding_limiter = WindingLimiter(
			graph.plane[p],
			graph.winding[p].clone(),
			None);
    	visited = recursive_leaf_flow(
    		&graph,
			p,
			graph.leaf_into[p],
			&mightsee,
			winding_limiter,
			&portalflood,
			&mut portalvis,
			&progress
    	);
    	//println!(" vis = {:?}", portalvis[p]);
    	progress[p].store(true, Ordering::Release);
    }
}

#[cfg(test)]
mod tests {
	use super::*;
    #[test]
    fn elim1() {
    	let mut deriv: BitfieldTracker = BitfieldTracker (0);
    	let mut ineqs: Vec<(Row, Relation, <BitfieldTracker as DeriviationTracker>::T)> = Vec::new();
    	ineqs.push(([ 1.0,  0.0,  0.0, 0.0, 0.0, 0.0], Relation::GeZero, deriv.new_item()));
    	ineqs.push(([-1.0, -2.0,  6.0, 0.0, 0.0, 0.0], Relation::GeZero, deriv.new_item()));
    	ineqs.push(([ 1.0,  1.0, -2.0, 0.0, 0.0, 0.0], Relation::GeZero, deriv.new_item()));
    	ineqs.push(([-1.0,  1.0,  3.0, 0.0, 0.0, 0.0], Relation::GeZero, deriv.new_item()));
    	ineqs.push(([ 0.0,  1.0,  0.0, 0.0, 0.0, 0.0], Relation::GeZero, deriv.new_item()));

    	let sol = any_solution(&ineqs, 0, 1, &mut deriv);
    	println!("sol {:?}", sol);

    	assert!(false);
    }

    #[test]
    fn elim2() {
    	let prt_test = "PRT1
3
2
4 0 1 (-1 -1 -1 ) (-1 -1 1 ) (-1 1 1 ) (-1 1 -1 )
4 1 2 (1 -1 -1 ) (1 -1 1 ) (1 1 1 ) (1 1 -1 )
";

	    let prt = prt_test;

    	let (nleafs, _, lines) = parse_prt(&prt).unwrap();
    	println!("lines {:?}", lines);

    	let graph = prtlines_to_graph(nleafs, &lines);
    	println!("graph {:?}", graph);

    	process_graph(&graph);

    	assert!(false);
    }
}
