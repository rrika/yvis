use crate::fm::DerivationTracker;

pub type N = f64;

#[derive(Clone, Debug)]
pub struct Winding { pub points: Vec<[N; 3]> }

#[derive(Copy, Clone, Debug, Default)]
pub struct Plane(pub N, pub N, pub N, pub N);

impl Plane {
	pub fn with_dist(&self, dist: N) -> Plane {
		Plane(self.0, self.1, self.2, dist)
	}
}

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

use crate::fm;

impl Winding {
	pub fn ray_passes(&self, _start: (N, N, N), _dir: (N, N, N)) -> bool {
		panic!("todo")
	}
}

pub fn fm_see_through_portals(portals: &Vec<&Winding>) -> bool
{
	let mut _portalid = 0;

	let num_constraints: usize = portals.iter().map(|p| p.points.len()).sum();

	let vis = if num_constraints >= 52 {
		let mut system: fm::System<Vec<u64>> = Vec::new();
	   	let mut deriv: fm::LongBitfieldTracker = fm::LongBitfieldTracker (0);
		for portal in portals {
			for i in 0..portal.points.len() {
				let b = portal.points[i];
				let a = portal.points[(i+1) % portal.points.len()];
				let x = ncross(a, b);
				let d = nsub(a, b);
				let row = [x[0], x[1], x[2], d[0], d[1], d[2]];
				system.push((row, fm::Relation::GtZero, deriv.new_item()));
			}
			_portalid += 1;
		}

		print!("[{:?} portals {:?} constraints]", portals.len(), system.len());
		let sol = fm::any_solution(&system, 0, 5, &mut deriv);
		if let Ok(_) = sol { true } else { false }
	} else {
		let mut system: fm::System<u64> = Vec::new();
	   	let mut deriv: fm::BitfieldTracker = fm::BitfieldTracker (0);
		for portal in portals {
			for i in 0..portal.points.len() {
				let b = portal.points[i];
				let a = portal.points[(i+1) % portal.points.len()];
				let x = ncross(a, b);
				let d = nsub(a, b);
				let row = [x[0], x[1], x[2], d[0], d[1], d[2]];
				system.push((row, fm::Relation::GtZero, deriv.new_item()));
			}
			_portalid += 1;
		}

		print!("[{:?} portals {:?} constraints]", portals.len(), system.len());
		let sol = fm::any_solution(&system, 0, 5, &mut deriv);
		if let Ok(_) = sol { true } else { false }
	};
	println!(" {:}", if vis { "vis" } else { "no" });
	vis
}

use lpsolve;
#[allow(unused_imports)]
use lpsolve_sys;

pub fn lpsolve_see_through_portals(portals: &Vec<&Winding>) -> bool {
	let num_constraints: usize = portals.iter().map(|p| p.points.len()).sum();
	let mut p = lpsolve::Problem::new(num_constraints as i32, 7).unwrap();
	unsafe { lpsolve_sys::set_verbose(p.to_lprec(), 0); } // shuddup

	let score: [f64; 8] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0];
	p.set_objective_function(&score[..]);
	p.set_unbounded(1);
	p.set_unbounded(2);
	p.set_unbounded(3);
	p.set_unbounded(4);
	p.set_unbounded(5);
	p.set_unbounded(6);

	for portal in portals {
		for i in 0..portal.points.len() {
			let b = portal.points[i];
			let a = portal.points[(i+1) % portal.points.len()];
			let x = ncross(a, b);
			let d = nsub(a, b);
			let row = [0.0, x[0], x[1], x[2], d[0], d[1], d[2], 1.0];
			//println!("row {:?}", row);
			p.add_constraint(&row[0..8], 0.5, lpsolve::ConstraintType::Ge);
		}
	}
	// println!("current size {:?}x{:?}", p.num_cols(), p.num_rows());
	// for i in 1..=p.num_rows() {
	// 	let mut row: [f64; 7] = [0.0; 7];
	// 	p.get_row(&mut row, i);
	// 	println!("rrw {:?} {:?}", row, p.get_constraint_type(i));
	// }

	let vis = match p.solve() {
		lpsolve::SolveStatus::Optimal => {
			let mut sol: [f64; 7] = [0.0; 7];
			p.get_solution_variables(&mut sol);
			//println!("lpsol {:?}", sol);
			sol[6] == 0.0
		},
		lpsolve::SolveStatus::Infeasible => false,
		lpsolve::SolveStatus::Unbounded => true, // ???
		status => {
			p.write_lp(&mut std::io::stdout());
			println!("{:?}", status);
			panic!()
		}
	};
	//println!(" {:}", if vis { "vis" } else { "no" });
	vis
}
