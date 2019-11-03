use crate::geometry::*;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Cuboid(N, N, N, N, N, N);
pub struct Rectangle(N, N, N, N);
struct Line(N, N, N);

impl Cuboid {
	pub fn project_x(&self) -> Rectangle { Rectangle(self.1, self.2, self.4, self.5) }
	pub fn project_y(&self) -> Rectangle { Rectangle(self.2, self.0, self.5, self.3) }
	pub fn project_z(&self) -> Rectangle { Rectangle(self.0, self.1, self.3, self.4) }
}

impl Line {
	pub fn through((a, b): (N, N), (c, d): (N, N)) -> Line {
		let e = d-b;
		let f = a-c;
		Line(e, f, -(a*e+b*f))
	}
	pub fn extend_x(&self) -> Plane { Plane(0.0, self.0, self.1, self.2) }
	pub fn extend_y(&self) -> Plane { Plane(self.1, 0.0, self.0, self.2) }
	pub fn extend_z(&self) -> Plane { Plane(self.0, self.1, 0.0, self.2) }
}

use float_ord::FloatOrd;
use crate::winding::Winding;

impl From<&Winding> for Cuboid {
    fn from(w: &Winding) -> Cuboid {
    	let xs = || w.points.iter().map(|[x, _, _]| FloatOrd(*x));
    	let ys = || w.points.iter().map(|[_, y, _]| FloatOrd(*y));
    	let zs = || w.points.iter().map(|[_, _, z]| FloatOrd(*z));
    	Cuboid(
    		xs().min().unwrap().0,
    		ys().min().unwrap().0,
    		zs().min().unwrap().0,
    		xs().max().unwrap().0,
    		ys().max().unwrap().0,
    		zs().max().unwrap().0
    	)
    }
}

impl PlaneSeparable<Cuboid> for Cuboid {
	fn separate<X, F: FnMut(Plane)->Option<X>>(&self, other: &Cuboid, mut f: F) -> Option<()> {
		if let Some((xa, xb)) = separate_rect(self.project_x(), other.project_x()) {
			f(xa.extend_x())?;
			f(xb.extend_x())?;
		}
		if let Some((ya, yb)) = separate_rect(self.project_y(), other.project_y()) {
			f(ya.extend_y())?;
			f(yb.extend_y())?;
		}
		if let Some((za, zb)) = separate_rect(self.project_z(), other.project_z()) {
			f(za.extend_z())?;
			f(zb.extend_z())?;
		}
		Some(())
	}
}

fn separate_rect(a: Rectangle, b: Rectangle) -> Option<(Line, Line)> {
	let classify = |i, j, k, l| {
		if j <= k { return 4 } // ijkl []()
		if l <= i { return 5 } // klij ()[]
		return 0
	};

	// 01 21
	// 03 23

	let (ia, ib, ic, id) = match (
		classify(a.0, a.2, b.0, b.2),
		classify(a.1, a.3, b.1, b.2))
	{
		(0, 0) => return None,
		(0, 4) => ((a.0, a.3), (b.2, b.1), (b.0, b.1), (a.2, a.3)),
		(0, 5) => ((a.2, a.1), (b.0, b.3), (b.2, b.3), (a.0, a.1)),

		(4, 0) => ((a.2, a.1), (b.0, b.3), (b.0, b.1), (a.2, a.3)),
		(4, 4) => ((a.0, a.3), (b.2, b.1), (b.0, b.3), (a.2, a.1)),
		(4, 5) => ((a.2, a.3), (b.0, b.1), (b.2, b.3), (a.0, a.1)),

		(5, 0) => ((a.0, a.1), (b.2, b.3), (b.2, b.1), (a.0, a.3)),
		(5, 4) => ((a.0, a.1), (b.2, b.3), (b.0, b.1), (a.2, a.3)),
		(5, 5) => ((a.2, a.1), (b.0, b.3), (b.2, b.1), (a.0, a.3)),

		(_, _) => panic!("invalid")
	};
	Some((Line::through(ib, ia), Line::through(id, ic)))
}

fn fmin(a: f64, b: f64) -> f64 {
	if a < b {a} else {b}
}

fn fmax(a: f64, b: f64) -> f64 {
	if a < b {b} else {a}
}

fn cut_four(Plane(x, y, z, w): Plane, a: N, b: N, d: N, e: N) -> N {
	let numbers = [
		a*x + b*y,
		d*x + b*y,
		a*x + e*y,
		d*x + e*y
	];
	let m = fmax(
		fmax(numbers[0], numbers[1]),
		fmax(numbers[2], numbers[3]));
	return (w-m)/z;
}

pub fn cut_cuboid(Plane(x, y, z, w): Plane, Cuboid(mut a, mut b, mut c, mut d, mut e, mut f): Cuboid) -> Cuboid {
	if z != 0.0 {
		let nz = cut_four(Plane(x, y, z, w), a, b, d, e);
		if z > 0.0 {
			c = fmax(c, nz);
		} else {
			f = fmin(f, nz)
		}
	}

	if x != 0.0 {
		let nz = cut_four(Plane(y, z, x, w), b, c, e, f);
		if x > 0.0 {
			a = fmax(a, nz);
		} else {
			d = fmin(d, nz)
		}
	}

	if y != 0.0 {
		let nz = cut_four(Plane(z, x, y, w), c, a, f, d);
		if y > 0.0 {
			b = fmax(b, nz);
		} else {
			e = fmin(e, nz)
		}
	}
	Cuboid(a, b, c, d, e, f)
}

impl PlaneChoppable for Cuboid {
	type Output = Cuboid;
	fn borrow(owned: &Cuboid) -> &Cuboid {
		&owned
	}
	fn chop(&self, plane: Plane) -> Option<Option<Cuboid>> {
		let nc = cut_cuboid(plane, *self);
		if nc.3 < nc.0 ||
		   nc.4 < nc.1 ||
		   nc.5 < nc.2
		{
			return None
		}
		Some(Some(nc))
	}
}

#[cfg(test)]
mod tests {
	use super::*;
    #[test]
    fn axis_chops() {
    	let c = Cuboid(-2.0, -2.0, -2.0, 2.0, 2.0, 2.0);

    	assert_eq!(c.chop(Plane( 1.0,  0.0,  0.0,  1.0)), Some(Some(Cuboid( 1.0, -2.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  1.0,  0.0,  1.0)), Some(Some(Cuboid(-2.0,  1.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  0.0,  1.0,  1.0)), Some(Some(Cuboid(-2.0, -2.0,  1.0,  2.0,  2.0,  2.0))));

    	assert_eq!(c.chop(Plane( 1.0,  0.0,  0.0, -1.0)), Some(Some(Cuboid(-1.0, -2.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  1.0,  0.0, -1.0)), Some(Some(Cuboid(-2.0, -1.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  0.0,  1.0, -1.0)), Some(Some(Cuboid(-2.0, -2.0, -1.0,  2.0,  2.0,  2.0))));

    	assert_eq!(c.chop(Plane(-1.0,  0.0,  0.0, -1.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  1.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0, -1.0,  0.0, -1.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  2.0,  1.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  0.0, -1.0, -1.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  2.0,  2.0,  1.0))));

    	assert_eq!(c.chop(Plane(-1.0,  0.0,  0.0,  1.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0, -1.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0, -1.0,  0.0,  1.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  2.0, -1.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  0.0, -1.0,  1.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  2.0,  2.0, -1.0))));

    	assert_eq!(c.chop(Plane( 2.0,  0.0,  0.0,  2.0)), Some(Some(Cuboid( 1.0, -2.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  2.0,  0.0,  2.0)), Some(Some(Cuboid(-2.0,  1.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 0.0,  0.0,  2.0,  2.0)), Some(Some(Cuboid(-2.0, -2.0,  1.0,  2.0,  2.0,  2.0))));
    }

    #[test]
    fn diagnoal_chops() {
    	let c = Cuboid(-2.0, -2.0, -2.0, 2.0, 2.0, 2.0);

    	assert_eq!(c.chop(Plane( 2.0,  1.0,  0.0,  2.0)), Some(Some(Cuboid( 0.0, -2.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane(-2.0, -1.0,  0.0, -2.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane( 2.0,  1.0,  0.0, -2.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  2.0,  2.0,  2.0))));
    	assert_eq!(c.chop(Plane(-2.0, -1.0,  0.0,  2.0)), Some(Some(Cuboid(-2.0, -2.0, -2.0,  0.0,  2.0,  2.0))));
    }
}
