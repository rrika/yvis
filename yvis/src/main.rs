use std::env;
use std::fs::File;
use std::io::Read;
use std::io::Write;

mod bsp;
mod prt;
mod pvsrle;
use prt::*;
use pvsrle::*;

fn write_ppm(name: &String, layers: &[([u8; 3], &Vec<Vec<bool>>)]) {
	let mut blob: Vec<u8> = Vec::new();
	let first_layer = layers[0].1;
	let h = first_layer.len();
	let w = first_layer[0].len();
	let header = format!("P6\n{} {}\n255\n", w, h);
	blob.extend_from_slice(header.as_ref());

	for y in 0..h {
		for x in 0..w {
			let mut color: [u8; 3] = [0; 3];
			for layer in layers {
				if layer.1[y][x] {
					color = layer.0;
					break;
				}
			}
			blob.extend_from_slice(&color[..]);
		}
	}

	let mut mfile = File::create(name).unwrap();
	mfile.write_all(&blob).expect("write failed");
}

fn load_prt(fname: &str) -> (usize, Vec<PRTLine>) {
	if fname.ends_with(".prt") {
		let (nleafs, _, lines) = {
			let mut file = File::open(fname).unwrap();
			let mut prt = String::new();
			file.read_to_string(&mut prt).unwrap();
			parse_prt(&prt).unwrap()
		};
		(nleafs, lines)
	} else if fname.ends_with(".bsp") {
		let bspdata = {
			let mut file = File::open(fname).unwrap();
			let mut blob = Vec::<u8>::new();
			file.read_to_end(&mut blob).unwrap();
			bsp::extract_relevant_data(&blob)
		};
		let lines = bsp::redo_prt(&bspdata);
		(bspdata.leafs.iter().map(|l|l.cluster as usize).filter(|c|*c!=!0).max().unwrap_or(0)+1, lines)
	} else {
		panic!("either a .prt or .bsp file is required")
	}
}

fn main() {
	let (fname, oname, mname) = {
		let mut args = env::args();
		args.next().unwrap();
		(args.next().expect("please provide a prt file"),
		 args.next().expect("please provide a pvs file"),
		 args.next())
	};

	let (nleafs, lines) = load_prt(&fname);

	//union_find_leafs(nleafs, &lines);
	let graph = prtlines_to_graph(nleafs, &lines);
	//let (pvs_fast, pvs) = process_graph(&graph);
	let (pvs_fast, pvs) = portalvis::report::report(&graph);

	let leaf_pvs_fast = graph.portalvis_to_leafvis(&pvs_fast);
	let leaf_pvs = graph.portalvis_to_leafvis(&pvs);

	if let Some(mname) = mname {
		write_ppm(&mname, &[
			([255, 255, 255], &leaf_pvs),
			//([80, 60, 60], &leaf_pvs_fast)
		]);
	}

	let cvis = compress_vis(&leaf_pvs, &leaf_pvs);

	let mut ofile = File::create(oname).unwrap();
	ofile.write_all(&cvis).expect("write failed");
}
