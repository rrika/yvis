use std::fs::*;
use std::io::Read;
use portalvis::*;

mod prt;
use prt::*;

fn main() {
	//let mut file = File::open("data/de_solace.prt").unwrap();
	let mut file = File::open("data/de_cinder.prt").unwrap();
	//let mut file = File::open("data/test.prt").unwrap();
	let mut prt_solace = String::new();
	file.read_to_string(&mut prt_solace).unwrap();
   	let (nleafs, _, lines) = parse_prt(&prt_solace).unwrap();
	union_find_leafs(nleafs, &lines);
   	let graph = prtlines_to_graph(nleafs, &lines);
   	//println!("{:?}", graph);
   	process_graph(&graph);
}
