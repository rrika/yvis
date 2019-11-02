use portalvis::graph::Portal;
use portalvis::geometry::Winding;
use portalvis::geometry::plane_through;
use portalvis::graph::LeafGraph;
use std::collections::HashMap;

#[derive(Clone, Debug, Default)]
pub struct PRTLine(usize, usize, Vec<[f32; 3]>);

pub fn parse_prt(data: &str) -> Result<(usize, usize, Vec<PRTLine>), &'static str> {

	let mut prtlines = Vec::<PRTLine>::new();

	let mut line_iter = data.split("\r\n");
	let magic = line_iter.next().ok_or("prt line 1 empty")?;
	let nleafs = line_iter.next().ok_or("prt line 2 empty")?;
	let nportals = line_iter.next().ok_or("prt line 3 empty")?;

	if magic != "PRT1" { return Err("first line isn't PRT1") }
	let nleafs_i = nleafs.parse::<usize>().map_err(|_| "nleafs not an integer")?;
	let nportals_i = nportals.parse::<usize>().map_err(|_| "nportals not an integer")?;

	let mut nleafs_implicit = 0;

	for _p in 0..nportals_i {
		let line = line_iter.next().ok_or("file truncated")?;
		let mut parts = line.split("(");
		let mut nlab = parts.next().ok_or("prt line without points")?.split(" ");
		let npoints = nlab.next().ok_or("prt line missing npoints")?;
		let leaf_a = nlab.next().ok_or("prt line missing leaf_a")?;
		let leaf_b = nlab.next().ok_or("prt line missing leaf_b")?;

		let mut winding = Vec::<[f32; 3]>::new();

		for _i in 0..npoints.parse::<usize>().map_err(|_| "npoints not an integer")? {
			// expected format is (x y z ) with space after the z
			let mut xyz = parts.next().ok_or("prt line with fewer points than declared")?.split(" ");
			let x = xyz.next().ok_or("prt point lacking x")?;
			let y = xyz.next().ok_or("prt point lacking y")?;
			let z = xyz.next().ok_or("prt point lacking z")?;

			winding.push([
				x.parse().map_err(|_| "prt point x not an integer")?,
				y.parse().map_err(|_| "prt point y not an integer")?,
				z.parse().map_err(|_| "prt point z not an integer")?
			]);
		}

		let leaf_a_i = leaf_a.parse::<usize>().map_err(|_| "prt leaf_a not an integer")?;
		let leaf_b_i = leaf_b.parse::<usize>().map_err(|_| "prt leaf_b not an integer")?;

		prtlines.push(PRTLine(leaf_a_i, leaf_b_i, winding));

		nleafs_implicit = std::cmp::max(nleafs_implicit, leaf_a_i+1);
		nleafs_implicit = std::cmp::max(nleafs_implicit, leaf_b_i+1);
	}

	Ok((nleafs_i, nleafs_implicit, prtlines))
}

pub fn union_find_leafs(nleafs: usize, lines: &Vec<PRTLine>) {
	let mut repr: Vec<usize> = Vec::new();
	for i in 0..nleafs {
		repr.push(i)
	}
	for PRTLine(src, dst, _) in lines {
		let mut rsrc = *src;
		let mut rdst = *dst;
		while rsrc != repr[rsrc] {rsrc = repr[rsrc]}
		while rdst != repr[rdst] {rdst = repr[rdst]}
		let r = std::cmp::min(rsrc, rdst);
		repr[*src] = r;
		repr[*dst] = r;
	}
	let mut m: HashMap<usize, usize> = HashMap::new();
	for rn in &repr {
		let mut n = *rn;
		while n != repr[n] {n = repr[n]}
		let counter = m.entry(n).or_insert(0);
		*counter += 1;
	}
	println!("nleafs: {:?}", nleafs);
	for (r, count) in m {
		println!("{:?}: {:?}", r, count);
	}
}

type N = f64;

pub fn prtlines_to_graph(nleafs: usize, lines: &Vec<PRTLine>) -> LeafGraph {
	let mut graph = LeafGraph::default();

	let nportals = 2 * lines.len(); // prt portals are two-sided
	graph.leaf_from.reserve(nleafs);
	graph.portals.reserve(nportals);
	for _i in 0..nleafs { graph.leaf_from.push(Vec::new()) }

	for PRTLine(src, dst, points) in lines {
		let mut points: Vec<[N; 3]> = points.iter().map(|a|[a[0] as N, a[1] as N, a[2] as N]).collect();
		let pi = graph.portals.len();
		graph.leaf_from[*src].push(pi+0);
		graph.leaf_from[*dst].push(pi+1);
		// println!("leaf[{}] += {}", src, pi+0);
		// println!("leaf[{}] += {}", dst, pi+1);
		points.reverse();
		let plane = plane_through(points[0], points[1], points[2]);
		let portal_fwd = Portal {
			winding: Winding{points: points.clone()},
			plane:   plane,
			leaf_from: *src,
			leaf_into: *dst
		};
		points.reverse();
		let portal_bwd = Portal {
			winding: Winding{points: points},
			plane:   -plane,
			leaf_from: *dst,
			leaf_into: *src
		};
		// println!("portals[{}] = {{from: {}, into {}}}", graph.portals.len(), portal_fwd.leaf_from, portal_fwd.leaf_into);
		graph.portals.push(portal_fwd);
		// println!("portals[{}] = {{from: {}, into {}}}", graph.portals.len(), portal_bwd.leaf_from, portal_bwd.leaf_into);
		graph.portals.push(portal_bwd);
	}

	graph
}

#[cfg(test)]
use portalvis::graph::process_graph;

#[test]
fn prttest() {
	let prt_test = "PRT1\r
4\r
3\r
4 0 1 (-1 1 -1 ) (-1 1 1 ) (-1 -1 1 ) (-1 -1 -1 )\r
4 1 2 (0 -0.5 -1 ) (0 -0.5 0.5 ) (0 -1 0.5 ) (0 -1 -1 )\r
4 2 3 (1 1 -1 ) (1 1 1 ) (1 -1 1 ) (1 -1 -1 )\r
";

    let prt = prt_test;

	let (nleafs, _, lines) = parse_prt(&prt).unwrap();
	println!("lines {:?}", lines);

	let graph = prtlines_to_graph(nleafs, &lines);
	println!("graph {:?}", graph);

	process_graph(&graph);

	println!("fm see 0");
	portalvis::geometry::lpsolve_see_through_portals(&vec![
			&graph.portals[0].winding
	]);
	println!("fm see 0-1");
	portalvis::geometry::lpsolve_see_through_portals(&vec![
			&graph.portals[0].winding,
			&graph.portals[1].winding
	]);
	println!("fm see 0-2");
	portalvis::geometry::lpsolve_see_through_portals(&vec![
			&graph.portals[0].winding,
			&graph.portals[2].winding,
	]);
	println!("fm see 2-4");
	portalvis::geometry::lpsolve_see_through_portals(&vec![
			&graph.portals[0].winding,
			&graph.portals[2].winding,
	]);
	println!("fm see 0-2-4");
	portalvis::geometry::lpsolve_see_through_portals(&vec![
			&graph.portals[0].winding,
			&graph.portals[2].winding,
			&graph.portals[4].winding
	]);

	assert!(false);
}
