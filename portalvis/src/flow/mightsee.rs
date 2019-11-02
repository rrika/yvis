use crate::bits::*;

pub enum MightseeResult<T> {
	Unchanged,
	Reduced(T),
	None
}

pub trait Mightsee<Filter, Confirm> : BitVecLike + Sized + Clone + Default where
	Filter: BitVecLike,
	Confirm: BitVecLike
{
	fn filter_for_unexplored<'a>(&'a self, f: &Filter, c: &Confirm) -> MightseeResult<Self>;
}

impl Mightsee<DenseBits, DenseBits> for DenseBits {
	fn filter_for_unexplored<'a>(&'a self, f: &DenseBits, c: &DenseBits) -> MightseeResult<DenseBits> {
		//let p_self = bitvec_pop(self);
		//let p_f = bitvec_pop(f);
		//let p_c = bitvec_pop(c);
		let (nmc, any_unexplored) = bitvec_and_and_not(self, f, c);
		//let p_nmc = bitvec_pop(&nmc);
		//println!("&&~ {} {} {} {}", p_self, p_f, p_c, p_nmc);
		if any_unexplored {
			MightseeResult::Reduced(nmc)
		} else {
			MightseeResult::None
		}
	}
}

impl<F, C> Mightsee<F, C> for SparseBits where
	F: BitVecLike,
	C: BitVecLike
{
	fn filter_for_unexplored<'a>(&'a self, f: &F, c: &C) -> MightseeResult<SparseBits> {
		//println!("&&~ sparse {}", values.len());
		let newvalues: Vec<usize>
			= self.0.iter().map(|v|*v).filter(|i: &usize| f[*i] && !c[*i]).collect();
		if newvalues.len() == 0 {
			MightseeResult::None
		} else if self.0.len() == newvalues.len() {
			MightseeResult::Unchanged
		} else {
			MightseeResult::Reduced(SparseBits(newvalues))
		}
	}
}
impl Mightsee<DenseBits, DenseBits> for ChunkedBits {
	fn filter_for_unexplored<'a>(&'a self, f: &DenseBits, c: &DenseBits) -> MightseeResult<ChunkedBits> {
		//println!("&&~ chunked {}", self.0.len());
		let newchunks: Vec<(usize, u64)>
			= self.0.iter().map(|(i, w)| (*i, *w & f.0[*i] & !c.0[*i])).filter(|(_, w)|*w != 0).collect();
		if  newchunks.len() > 0 {
			MightseeResult::Reduced(ChunkedBits(newchunks))
		} else {
			MightseeResult::None
		}
	}
}
