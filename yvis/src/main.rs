use std::env;
use std::fs::File;
use std::io::Read;
use std::io::Write;
use portalvis::*;

mod prt;
mod pvsrle;
use prt::*;
use pvsrle::*;

fn main() {
	let (fname, oname) = {
		let mut args = env::args();
		args.next().unwrap();
		(args.next().expect("please provide a prtfile"),
		 args.next())
	};

	let (nleafs, _, lines) = {
		let mut file = File::open(fname).unwrap();
		let mut prt = String::new();
		file.read_to_string(&mut prt).unwrap();
	   	parse_prt(&prt).unwrap()
	};

	union_find_leafs(nleafs, &lines);
   	let graph = prtlines_to_graph(nleafs, &lines);
   	let pvs = process_graph(&graph);

   	let cvis = compress_vis(&pvs, &pvs);

   	if let Some(roname) = oname {
   		let mut ofile = File::create(roname).unwrap();
   		ofile.write_all(&cvis).expect("write failed");
   	}
}
