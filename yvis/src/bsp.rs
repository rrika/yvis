use portalvis::winding::chop_winding;
use std::collections::HashMap;
use portalvis::geometry::nsub;
use portalvis::geometry::nadd;
use portalvis::geometry::N;
use portalvis::geometry::nmul;
use portalvis::geometry::ncross;
use portalvis::winding::Winding;
use portalvis::geometry::Plane;
use std::io::Cursor;
use std::io::Seek;
use byteorder::{LittleEndian, ReadBytesExt};

// dleaf_t

pub struct BspNode { plane: i32, child0: i32, child1: i32 }
pub struct BspPlane(f32, f32, f32, f32);
pub struct BspLeaf { contents: i32, cluster: i16 }
pub struct BspData {
	pub nodes:  Vec<BspNode>,
	pub planes: Vec<BspPlane>,
	pub leafs:  Vec<BspLeaf>
}

fn extract_node(reader: &mut Cursor<&[u8]>) -> Option<BspNode> {
	let plane  = (*reader).read_i32::<LittleEndian>().ok()?;
	let child0 = reader.read_i32::<LittleEndian>().unwrap();
	let child1 = reader.read_i32::<LittleEndian>().unwrap();
	reader.seek(std::io::SeekFrom::Current(20)).unwrap();
	Some(BspNode { plane, child0, child1 })
}

fn extract_plane(reader: &mut Cursor<&[u8]>) -> Option<BspPlane> {
	let x = reader.read_f32::<LittleEndian>().ok()?;
	let y = reader.read_f32::<LittleEndian>().unwrap();
	let z = reader.read_f32::<LittleEndian>().unwrap();
	let d = reader.read_f32::<LittleEndian>().unwrap();
	reader.seek(std::io::SeekFrom::Current(4)).unwrap();
	Some(BspPlane(x, y, z, d))
}

fn extract_leaf(reader: &mut Cursor<&[u8]>) -> Option<BspLeaf> {
	let contents = reader.read_i32::<LittleEndian>().ok()?;
	let cluster  = reader.read_i16::<LittleEndian>().ok()?;
	reader.seek(std::io::SeekFrom::Current(26)).unwrap();
	Some(BspLeaf { contents, cluster })
}

pub fn extract_relevant_data(bsp: &[u8]) -> BspData {
	let get_leaf = |index: usize| {
		let header_offset = 8 + 16*index;
		let mut reader = Cursor::new(&bsp[header_offset..header_offset+16]);
		let offset = reader.read_u32::<LittleEndian>().unwrap() as usize;
		let length = reader.read_u32::<LittleEndian>().unwrap() as usize;
		&bsp[offset..offset+length]
	};
	//println!("{:?}", &bsp[..16]);
	let mut nodes:  Vec<BspNode> = Vec::new();
	let mut planes: Vec<BspPlane> = Vec::new();
	let mut leafs:  Vec<BspLeaf> = Vec::new();
	{
		let mut node_cursor = Cursor::new(get_leaf(5));
		while let Some(node) = extract_node(&mut node_cursor) {
			nodes.push(node)
		}
		let mut plane_cursor = Cursor::new(get_leaf(1));
		while let Some(plane) = extract_plane(&mut plane_cursor) {
			planes.push(plane)
		}
		let mut leaf_cursor = Cursor::new(get_leaf(10));
		while let Some(leaf) = extract_leaf(&mut leaf_cursor) {
			leafs.push(leaf)
		}
	}
	BspData { nodes, planes, leafs }
}

// use portalvis::graph::Portal;
// use portalvis::winding::Winding;

// pub struct Portal {
// 	pub winding: Winding,
// 	pub plane:   Plane,
// 	pub leaf_from: usize,
// 	pub leaf_into: usize
// }

fn winding_from_plane(plane: Plane, size: i32) -> Winding {
	let n = [plane.0, plane.1, plane.2];
	let mut a = ncross(n, [0.0, 1.0, 0.0]);
	if a == [0.0, 0.0, 0.0] {
		a = ncross(n, [1.0, 0.0, 0.0])
	}
	let mut b = ncross(n, a);
	let l = (plane.0*plane.0 + plane.1*plane.1 + plane.2*plane.2).sqrt();
	let p = nmul(n, plane.3/l);
	a = nmul(a, size as N); 
	b = nmul(b, size as N);
	let v1 = nsub(nadd(p, a), b);
	let v2 = nsub(nsub(p, a), b);
	let v3 = nadd(nsub(p, a), b);
	let v4 = nadd(nadd(p, a), b);
	Winding{points: vec![v2, v3, v4, v1]}
}

enum Split {
	Pos,
	Neg,
	Across(Winding, Winding)
}

fn double_chop_winding(w: &[[N; 3]], plane: Plane) -> Split {
	let pos = chop_winding(w,  plane);
	let neg = chop_winding(w, -plane);
	match (pos, neg) {
		(Some(None),    None)          => Split::Pos,
		(None,          Some(None))    => Split::Neg,
		(Some(Some(a)), Some(Some(b))) => Split::Across(a, b),

		// be a bit lenient
		(Some(None),    _)             => Split::Pos,
		(_,             Some(None))    => Split::Neg,
		(Some(Some(_)), None)          => Split::Pos,
		(None,          Some(Some(_))) => Split::Neg,

		(pos, neg) => {
			println!("{:?}", w);
			println!("{:?}", plane);
			println!("{:?}", pos);
			println!("{:?}", neg);
			panic!("oh no")
		}
	}
}

fn replace_in_adjacent(adj: &mut [(Winding, Plane, i32, u32)], matchpi: u32,
	new_neighbour: i32, new_winding_o: Option<Winding>)
{
	for (w, _pl, n, pi) in adj {
		if *pi == matchpi {
			*n = new_neighbour;
			if let Some(new_winding) = new_winding_o {
				*w = new_winding;
			}
			return
		}
	}
	panic!("none found");
}

fn split_list(
	plane: Plane,
	olist: &[(Winding, Plane, i32, u32)],
	child0: i32,
	child1: i32,
	node_portals: &mut HashMap<i32, Vec<(Winding, Plane, i32, u32)>>,
	pcounter: &mut u32,
	depth: usize) ->
		(Vec<(Winding, Plane, i32, u32)>, Vec<(Winding, Plane, i32, u32)>)
{
	let mut plist: Vec<(Winding, Plane, i32, u32)> = Vec::new();
	let mut nlist: Vec<(Winding, Plane, i32, u32)> = Vec::new();
	for (w, p, n, pi) in olist {
		// println!("{:indent$} double_chop neighbour={} portalid={}, normal={:?}, distances={:?}",
		// 	"", n, pi, p, w.points.iter().map(|point| *p * *point).collect::<Vec<_>>(),
		// 	indent=depth*2);
		//println!("{:indent$}double_chop {:?} {:?} {:?} {:?}", "", w, p, n, pi, indent=depth*2);
		match double_chop_winding(&w.points[..],  plane) {
			Split::Pos => {
				//println!("{:indent$}  pos", "", indent=depth*2);
				plist.push((w.clone(), *p, *n, *pi));
				replace_in_adjacent(node_portals.entry(*n).or_default(), *pi, child0, None)
			},
			Split::Neg => {
				//println!("{:indent$}  neg", "", indent=depth*2);
				nlist.push((w.clone(), *p, *n, *pi));
				replace_in_adjacent(node_portals.entry(*n).or_default(), *pi, child1, None)
			},
			Split::Across(mut a, mut b) => {
				let pi2 = *pcounter; *pcounter += 1;
				//println!("{:indent$}  across {}/{}", "", pi, pi2, indent=depth*2);
				plist.push((a.clone(), *p, *n, *pi));
				nlist.push((b.clone(), *p, *n, pi2));
				a.points.reverse();
				b.points.reverse();
				replace_in_adjacent(node_portals.entry(*n).or_default(), *pi, child0, Some(a));
				node_portals.entry(*n).or_default().push((b, -*p, child1, pi2))
			}
		}
	}
	(plist, nlist)
}

fn boundary_portals(bplane: Plane, list: &[(Winding, Plane, i32, u32)]) -> Option<(Winding, Winding)> {
	let mut pw = Some(winding_from_plane(bplane, 100000));
	for (_, p, _, _) in list {
		pw = if let Some(pws) = pw {
			match chop_winding(&pws.points[..], *p) {
				Some(Some(w)) => Some(w),
				Some(None) => Some(pws),
				None => None,
			}
		} else { None }
	}
	let pw = pw?;
	let mut nw = pw.clone();
	nw.points.reverse();
	Some((pw, nw))
}

fn trace_tree(bspdata: &BspData, node: i32, depth: usize,
	node_portals: &mut HashMap<i32, Vec<(Winding, Plane, i32, u32)>>,
	pcounter: &mut u32) {
	if node < 0 {
		// println!("{:indent$}leaf {:?} contents={}", "",
		// 	!node, bspdata.leafs[(!node) as usize].contents, indent=depth*2);
	} else {
		let n = &bspdata.nodes[node as usize];
		let pl = {let BspPlane(x, y, z, d) = bspdata.planes[n.plane as usize]; Plane(x as N, y as N, z as N, d as N)};
		// println!("{:indent$}{:?} (plane[{}]={:?})", "", node, n.plane, pl, indent=depth*2);
		let olist = node_portals.remove(&node).unwrap();
		let (mut plist, mut nlist) = split_list(pl, &olist, n.child0, n.child1, node_portals, pcounter, depth);
		if let Some((pw, nw)) = boundary_portals(pl, &olist) {
			plist.push((pw,  pl, n.child1, *pcounter));
			nlist.push((nw, -pl, n.child0, *pcounter));
			*pcounter += 1;
		}
		node_portals.insert(n.child0, plist);
		node_portals.insert(n.child1, nlist);
		trace_tree(bspdata, n.child0, depth+1, node_portals, pcounter);
		trace_tree(bspdata, n.child1, depth+1, node_portals, pcounter);
	}
}

use crate::prt::PRTLine;

pub fn redo_prt(bspdata: &BspData) -> Vec<PRTLine> {
	let mut node_portals: HashMap<i32, Vec<(Winding, Plane, i32, u32)>> = HashMap::new();
	let mut pcounter: u32 = 0;
	node_portals.entry(0).or_default();
	trace_tree(bspdata, 0, 0, &mut node_portals, &mut pcounter);
	let mut prtlines: Vec<PRTLine> = Vec::new();
	for (k, v) in node_portals {
		for (w, _p, n, _i) in v {
			let mut k = (!k) as usize;
			let mut n = (!n) as usize;
			if bspdata.leafs[k].cluster == -1 { continue }
			if bspdata.leafs[n].cluster == -1 { continue }
			k = bspdata.leafs[k].cluster as usize;
			n = bspdata.leafs[n].cluster as usize;
			if k < n {
				prtlines.push(PRTLine(
					k,
					n,
					w.points))
			}
		}
	}
	prtlines
}
