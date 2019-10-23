fn bitpack_row(row: &Vec<bool>) -> Vec<u8> {
	let mut packed: Vec<u8> = Vec::new();
	let mut i = 0usize;
	let mut p = 0u8;
	while i < row.len() {
		p |= (row[i] as u8) << (i % 8);
		if i % 8 == 7 {
			packed.push(p);
			p = 0;
		}
		i += 1;
	}
	if i % 8 != 7 {
		packed.push(p)
	}
	packed
}

fn compress_row(dst: &mut Vec<u8>, row: &Vec<u8>) -> usize {
	let original_len = dst.len();
	let mut zeroes: u8 = 0;
	for i in row {
		if *i != 0 || zeroes == 255 {
			if zeroes > 0 { dst.push(zeroes) }
			dst.push(*i);
			zeroes = 0;
		}
		if *i == 0 {
			zeroes += 1;
		}
	};
	if zeroes > 0 { dst.push(zeroes) }
	original_len
}

pub fn compress_vis(pvs: &Vec<Vec<bool>>, pas: &Vec<Vec<bool>>) -> Vec<u8> {
	let n = pvs.len();
	let mut header: Vec<u32> = Vec::new(); header.resize(1 + 2*n, 0u32);
	let mut blob: Vec<u8> = Vec::new(); blob.resize(4 + 8*n, 0u8);
	header[0] = n as u32;
	for i in 0..n {
		let row = bitpack_row(&pvs[i]);
		header[2*i+1] = compress_row(&mut blob, &row) as u32;
	}
	for i in 0..n {
		let row = bitpack_row(&pas[i]);
		header[2*i+2] = compress_row(&mut blob, &row) as u32;
	}
	for i in 0..header.len() {
		blob[4*i + 0] = (header[i] >>  0) as u8;
		blob[4*i + 1] = (header[i] >>  8) as u8;
		blob[4*i + 2] = (header[i] >> 16) as u8;
		blob[4*i + 3] = (header[i] >> 24) as u8;
	}
	blob
}
