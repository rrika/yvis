use std::fs::*;
use std::io::Read;
use portalvis::*;

fn main() {
	let mut file = File::open("de_solace.prt").unwrap();
	//let mut file = File::open("test.prt").unwrap();
	let mut prt_solace = String::new();
	file.read_to_string(&mut prt_solace).unwrap();
   	let (nleafs, _, lines) = parse_prt(&prt_solace).unwrap();
   	let graph = prtlines_to_graph(nleafs, &lines);
   	process_graph(&graph);
}
