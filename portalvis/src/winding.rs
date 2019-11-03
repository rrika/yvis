
use crate::geometry::PlaneChoppable;
use crate::geometry::*;
use crate::fm;
use crate::fm::DerivationTracker;

#[derive(Clone, Debug)]
pub struct Winding{pub points: Vec<[N; 3]>}

impl Winding {
	pub fn ray_passes(&self, _start: (N, N, N), _dir: (N, N, N)) -> bool {
		panic!("todo")
	}
}

impl PlaneChoppable for [[N; 3]] {
	type Output = Winding;
	fn borrow(owned: &Winding) -> &[[N; 3]] {
		&owned.points
	}
	fn chop(&self, plane: Plane) -> Option<Option<Winding>> {
		chop_winding(self, plane)
	}
}

pub fn chop_winding(w: &[[N; 3]], p: Plane) -> Option<Option<Winding>> {
	//println!("chop_winding {:?} {:?}", p, w);
	let mut not_all = false;
	for point in w {
		if p * *point < 0.0 {
			not_all = true;
		}
	}
	if !not_all {
		return Some(None)
	}

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
	let mut i = w.len()-1;
	for j in 0..w.len() {
		let c = w[i];
		let d = w[j];
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
		Some(Some(Winding{points: nw}))
	} else {
		None
	}
}


pub fn clip_hourglass(source: &[[N; 3]], pass: &[[N; 3]], target: &[[N; 3]]) -> Option<Option<Winding>> {
	// None    None
	// All     Some(None)
	// Some(x) Some(Some(x))
	let mut target_clipped: Option<Winding> = None;

	for i in 0..source.len() {
		let j = (i+1) % source.len();
		for k in 0..pass.len() {
			let l = (k+1) % pass.len();

			let p1 = plane_through(
				source[i],
				source[j],
				pass[k]);

			let p2 = plane_through(
				source[i],
				pass[k],
				pass[l]);


			// all_below(&source, p1) is true
			if all_above(&pass, p1) {
				target_clipped = chop_winding(
					if let Some(tw) = &target_clipped { &tw.points[..] } else { target },
					p1)?.or(target_clipped);
			}

			// all_above(&pass.points, p2) is true
			if all_below(&source, p2) {
				target_clipped = chop_winding(
					if let Some(tw) = &target_clipped { &tw.points[..] } else { target },
					p2)?.or(target_clipped);
			}
		}
	}
	Some(target_clipped)
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
