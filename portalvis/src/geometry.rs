type N = f64;

#[derive(Clone, Debug)]
pub struct Winding { pub points: Vec<[N; 3]> }

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

pub fn fm_see_through_portals(portals: &Vec<&Winding>) -> bool
{
	let mut system: fm::System<u64> = Vec::new();
	let mut portalid = 0;

	for portal in portals {
		for i in 0..portal.points.len() {
			let b = portal.points[i];
			let a = portal.points[(i+1) % portal.points.len()];
			let x = ncross(a, b);
			let d = nsub(a, b);
			let row = [x[0], x[1], x[2], d[0], d[1], d[2]];
			system.push((row, fm::Relation::GtZero, 1u64 << portalid));
		}
		portalid += 1;
	}

   	let mut deriv: fm::BitfieldTracker = fm::BitfieldTracker (0);
	let sol = fm::any_solution(&system, 0, 5, &mut deriv);
	println!("see sol {:?}", sol);
	println!("");
	if let Ok(_) = sol { true } else { false }
}
