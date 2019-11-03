
use std::fmt::Debug;

pub type N = f64;

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

pub trait PlaneSeparable<T> {
	fn separate<X, F: FnMut(Plane)->Option<X>>(&self, other: &T, f: F) -> Option<()>;
}

pub trait PlaneChoppable {
	type Output;
	fn borrow(owned: &Self::Output) -> &Self;
	fn chop(&self, plane: Plane) -> Option<Option<Self::Output>>;
}

pub fn generalized_clip_hourglass<T: PlaneSeparable<T> + Debug, U: PlaneChoppable + Debug>(
	source: &T,
	pass:   &T,
	target: &U

) -> Option<Option<U::Output>> where U::Output: Debug {

	let mut target_clipped: Option<U::Output> = None;
	//println!("");
	{
		source.separate(pass, |plane| {
			//println!("plane: {:?}", plane);
			if let Some(x) = (if let Some(tw) = &target_clipped { U::borrow(tw) } else { target }).chop(plane)? {
				target_clipped = Some(x)
			}
			Some(())
		});
	}
	//println!("done");

	Some(target_clipped)
}
