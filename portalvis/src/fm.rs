use std::fmt::Debug;
use std::collections::HashMap;

type N = f64;
pub(crate) type Row = [N; 6];

#[derive(Copy, Clone, Debug)]
pub enum Relation {
    GeZero,
    GtZero
}

pub type System<T> = Vec<(Row, Relation, T)>;

#[derive(Debug)]
pub struct ReconstructionInfo<T: DeriviationTracker + Debug> {
    pub gt: Vec<(Row, T::T)>,
    pub ge: Vec<(Row, T::T)>,
    pub lt: Vec<(Row, T::T)>,
    pub le: Vec<(Row, T::T)>,
}

pub fn dot(lhs: Row, rhs: Row) -> N {
    lhs[0] * rhs[0] + 
    lhs[1] * rhs[1] + 
    lhs[2] * rhs[2] + 
    lhs[3] * rhs[3] + 
    lhs[4] * rhs[4] + 
    lhs[5] * rhs[5]
}

pub fn normalize_row(row: Row) -> Row {
    let len = dot(row, row).sqrt();
    if len == 0.0 { row } else { divide_row(row, len) }
}

pub fn reconstruct<T: DeriviationTracker+Debug>(lhs: Row, recon: &ReconstructionInfo<T>, _debugname: &str, deriv: &mut T)
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
    //println!("{:?} ... {} ... {:?}", l, debugname, g);
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
    fn either(&mut self, a: Self::T, b: Self::T) -> Self::T;
    fn empty(&mut self) -> Self::T;
    fn partition(&mut self, a: Self::T, b: Self::T) -> (Self::T, Self::T);
    fn count(&self, a: Self::T) -> u32;

    fn as_u64(&self, a: Self::T) -> u64;
}

#[derive(Debug)]
pub struct BitfieldTracker(pub(crate) usize);
impl DeriviationTracker for BitfieldTracker {
    type T = u64;
    fn new_item(&mut self) -> u64 {
        self.0 += 1;
        if self.0 >= 64 { panic!("out of bits") }
        1u64 << self.0-1
    }
    fn combine(&mut self, a: u64, b: u64) -> u64 { a | b }
    fn either(&mut self, a: u64, _b: u64) -> u64 { a } // could pick the one with fewer bits
    fn empty(&mut self) -> u64 { 0u64 }
    fn partition(&mut self, a: u64, b: u64) -> (u64, u64) { (a & b, a & !b) }
    fn count(&self, a: u64) -> u32 { a.count_ones() }

    fn as_u64(&self, a: u64) -> u64 { a }
}

pub fn fourier_motzkin_elimination<T: DeriviationTracker+Debug>(
    ineqs: &System<T::T>,
    axis: usize,
    deriv: &mut T,
    t_elim: T::T,
    t_axis: T::T,
    t_nonzeroes: &[T::T]
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

    let mut dedup = |rows: &Vec<(Row, T::T)>| -> Vec<(Row, T::T)> {
        let mut m: HashMap<[u64; 6], T::T> = HashMap::new();
        unsafe {
            for row in rows {
                m.entry(std::mem::transmute(normalize_row(row.0)))
                //m.entry(std::mem::transmute(row.0))
                    .and_modify(|t| *t = deriv.either(*t, row.1))
                    .or_insert(row.1);
            }
            m.iter().map(|(k, v)| (std::mem::transmute::<[u64; 6], Row>(*k), *v)).collect()
        }
    };

    reconstruct.gt = dedup(&reconstruct.gt);
    reconstruct.ge = dedup(&reconstruct.ge);
    reconstruct.lt = dedup(&reconstruct.lt);
    reconstruct.le = dedup(&reconstruct.le);

    let t_nonzeroes_mask = {
        let mut tc = deriv.empty();
        for t in t_nonzeroes { tc = deriv.combine(*t, tc) }
        tc
    };

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

                let t_ab = deriv.combine(*at, *bt);
                let t_abe = deriv.combine(t_ab, t_axis);

                // imbert acceleration
                let mut t_zero = deriv.empty();
                for axis in 0..6 {
                    if new_row[axis] == 0.into() {
                        t_zero = deriv.combine(t_zero, t_nonzeroes[axis])
                    }
                }
                t_zero = deriv.partition(t_zero, t_elim).0;
                let (e, hz) = deriv.partition(t_abe, t_elim);
                let (z, h) = deriv.partition(hz, t_nonzeroes_mask);
                let (io, _) = deriv.partition(z, t_zero);
                let lhs = 1 + deriv.count(e);
                let mid = deriv.count(h);
                if lhs != mid {
                    let ei = deriv.combine(e, io);
                    let rhs = 1 + deriv.count(ei);
                    //println!("{:?} {:?} {:?} E={:016b}  Z={:016b}  H={:016b}  a={:016b} b={:016b}", lhs, mid, rhs,
                    //    deriv.as_u64(e),
                    //    deriv.as_u64(z),
                    //    deriv.as_u64(h),
                    //    deriv.as_u64(*at),
                    //    deriv.as_u64(*bt));
                    if lhs > mid || mid > rhs {
                        continue
                    }
                }
                new_ineqs.push((new_row, *rel, t_abe));
            }
        }
    }
    (reconstruct, new_ineqs)
}

pub fn any_solution<T: DeriviationTracker+Debug>(
    ineqs: &System<T::T>,
    first_axis: usize,
    last_axis: usize,
    deriv: &mut T) -> Result<Row, T::T>
{
    let mut ineqs = ineqs.clone();
    let mut reconstruction: Vec<(usize, ReconstructionInfo<T>)> = Vec::new();

    let mut nonzeroes: Vec<T::T> = Vec::new();
    for axis in 0..6 {
        let nonzero = deriv.new_item();
        nonzeroes.push(nonzero);
        for (row, _relation, t) in &mut ineqs {
            if row[axis] != 0.into() {
                *t = deriv.combine(*t, nonzero);
            }
        }
    }

    let mut t_elim = deriv.empty();

    for axis in first_axis..=last_axis {
        let t_axis = deriv.new_item();
        t_elim = deriv.combine(t_elim, t_axis);

        //println!("axis {:?}, num_ineqs {:?}", axis, ineqs.len());
        //for ineq in &ineqs { println!("{:?}", ineq); } println!("--");
        reconstruction.push((axis, {
            let (recon, new_ineqs) = fourier_motzkin_elimination(
                &ineqs, axis, deriv,
                t_elim,
                t_axis,
                &nonzeroes[..]);
            ineqs = new_ineqs;
            recon
        }));
        print!("-{:?}", ineqs.len());
    }
    let mut solution: Row = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    for (axis, rinfo) in reconstruction.iter().rev() {
        solution[*axis] = reconstruct::<T>(solution, &rinfo, &["x", "y", "z", "a", "b", "c"][*axis], deriv)?;
        //println!("sol {:?}", solution);
    }
    Ok(solution)
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
    }
}
