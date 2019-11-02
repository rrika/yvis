use std::ops;

pub trait BitVecLike: Sized + ops::Index<usize, Output=bool> {}
impl<T> BitVecLike for T where T: Sized + ops::Index<usize, Output=bool> {}


#[derive(Clone, Default)]
pub struct SparseBits(pub Vec<usize>);

#[derive(Clone, Default)]
pub struct ChunkedBits(pub Vec<(usize, u64)>);

#[derive(Clone, Default)]
pub struct DenseBits(pub Vec<u64>);

const SOMETRUE: bool = true;
const SOMEFALSE: bool = false;

impl ops::Index<usize> for SparseBits {
	type Output = bool;
	fn index(&self, index: usize) -> &bool {
		if self.0.binary_search(&index).is_ok() { &SOMETRUE } else { &SOMEFALSE }
	}
}

impl ops::Index<usize> for ChunkedBits {
	type Output = bool;
	fn index(&self, index: usize) -> &bool {
		let chunkindex = index / 64;
		let chunk: u64 = if let Ok(index) = self.0.binary_search_by_key(&chunkindex, |&(a,_)| a) {
			self.0[index].1
		} else { 0 };
		if 1 == 1 & (chunk >> (index % 64)) { &SOMETRUE } else { &SOMEFALSE }
	}
}

impl ops::Index<usize> for DenseBits {
	type Output = bool;
	fn index(&self, index: usize) -> &bool {
		let chunkindex = index / 64;
		let chunk = self.0[chunkindex];
		if 1 == 1 & (chunk >> (index % 64)) { &SOMETRUE } else { &SOMEFALSE }
	}
}

impl ops::BitOr<&DenseBits> for &DenseBits {
	type Output = DenseBits;

	fn bitor(self, b: &DenseBits) -> DenseBits {
		assert_eq!(self.0.len(), b.0.len());
		let mut o = Vec::new();
		o.reserve(self.0.len());
		for i in 0..self.0.len() {
			let v = self.0[i] | b.0[i];
			o.push(v);
		}
		DenseBits(o)
	}
}

pub trait SetBit {
	fn set(&mut self, i: usize);
}

impl SetBit for DenseBits {
	fn set(&mut self, i: usize) {
		self.0[i / 64] |= 1 << (i % 64)
	}
}

impl DenseBits {
	pub fn clear(&mut self, i: usize) {
		self.0[i / 64] &= !(1 << (i % 64))
	}

	pub fn count_ones(&self) -> usize {
		self.0.iter().map(|v|u64::count_ones(*v)).sum::<u32>() as usize
	}
}

pub fn bitvec_and_not_pop(a: &DenseBits, b: &DenseBits) -> usize {
	let mut c = 0usize;
	for i in 0..a.0.len() {
		c += (a.0[i] & !b.0[i]).count_ones() as usize;
	}
	c
}

pub fn bitvec_and_and_not(a: &DenseBits, b: &DenseBits, c: &DenseBits) -> (DenseBits, bool) {
	assert_eq!(a.0.len(), b.0.len());
	assert_eq!(a.0.len(), c.0.len());
	let mut o = Vec::new();
	let mut any = false;
	o.reserve(a.0.len());
	for i in 0..a.0.len() {
		let v = a.0[i] & b.0[i] & !c.0[i];
		o.push(v);

		if v != 0 { any = true }
	}
	return (DenseBits(o), any)
}

pub fn sparse_convert(bitvec: &DenseBits, n: usize) -> SparseBits {
	let mut v = Vec::<usize>::new();
	for i in 0..n {
		if bitvec[i] { v.push(i) }
	}
	SparseBits(v)
}

pub fn chunked_convert(bitvec: &DenseBits) -> ChunkedBits {
	ChunkedBits(bitvec.0.iter().cloned().enumerate().filter(|(_, w)| *w != 0).collect())
}
