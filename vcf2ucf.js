function b8_parse_vcf2(t) // t = vcf_line.split("\t")
{
	if (t.length < 6) return null;
	var a = [], pos = parseInt(t[1]) - 1;
	t[3] = t[3].toUpperCase(); t[4] = t[4].toUpperCase();
	var s = t[4].split(","); // list of ALT alleles
	// get CIGAR for freebayes
	var m3 = /CIGAR=([^;\t]+)/.exec(t[7]);
	var cigar = m3 != null? m3[1].split(",") : [];
	if (cigar.length && cigar.length != s.length) throw Error("Inconsistent ALT and CIGAR");
	// loop through each ALT allele
	for (var i = 0; i < s.length; ++i) {
		if (t[3].length == 1 && s[i].length == 1) { // SNP
			if (t[3] != s[i]) a.push([pos, pos+1, 0, t[3], s[i], i]);
		} else if (cigar.length) { // MNP or INDEL from freebayes
			var x = 0, y = 0;
			var m4, re = /(\d+)([MXID])/g;
			while ((m4 = re.exec(cigar[i])) != null) {
				var l = parseInt(m4[1]);
				if (m4[2] == 'X') {
					for (var j = 0; j < l; ++j) {
						var u = t[3].substr(x+j, 1), v = s[i].substr(y+j, 1);
						a.push([pos + x, pos+x+1, 0, u, v, i]);
					}
					x += l, y += l;
				} else if (m4[2] == 'I') {
					if (x == 0 || y == 0) throw Error("Leading I/D");
					var u = t[3].substr(x-1, 1), v = s[i].substr(y-1, l+1);
					a.push([pos + x - 1, pos+x, l, u, v, i]);
					y += l;
				} else if (m4[2] == 'D') {
					if (x == 0 || y == 0) throw Error("Leading I/D");
					var u = t[3].substr(x-1, l+1), v = s[i].substr(y-1, 1);
					a.push([pos + x - 1, pos+x+l, -l, u, v, i]);
					x += l;
				} else if (m4[2] == 'M') x += l, y += l;
			}
		} else { // MNP or INDEL from Platypus and others
			var d = s[i].length - t[3].length, rs, as;
			if (d > 0) {
				a.push([pos, pos + 1, d, t[3].charAt(0), s[i].substr(0, d + 1), i]);
				rs = 1, as = d + 1;
			} else if (d < 0) {
				a.push([pos, pos + 1 + (t[3].length + d), d, t[3].substr(0, -d + 1), s[i].charAt(0), i]);
				rs = -d + 1, as = 1;
			} else rs = as = 0;
			for (var j = 0; j < t[3].length - rs; ++j) { // check the rest
				var u = t[3].substr(rs + j, 1), v = s[i].substr(as + j, 1);
				if (u != v) a.push([pos + j + rs, pos + j + rs + 1, 0, u, v, i]);
			}
		}
	}
	return a; // [start, end, indelLen, ref, alt, i]
}

function b8_vcf2ucf(args)
{
	if (args.length == 0) {
		print("Usage: k8 vcf2ucf.js <in.vcf>");
		exit(0);
	}

	var file = args[0] == '-'? new File() : new File(args[0]);
	var buf = new Bytes();
	var srt_buf = [];

	while (file.readline(buf) >= 0) {
		if (buf.length == 0) continue; // skip empty lines
		var line = buf.toString();
		if (line.charAt(0) == '#') { // VCF header
			if (line.length > 1 && line.charAt(1) != '#')
				print('##ALT=<ID=M,Description="miscellaneous allele">');
			print(line);
		} else {
			var t = line.split("\t");
			t[1] = parseInt(t[1]);
			// print srt_buf
			while (srt_buf.length && (srt_buf[0][0] != t[0] || srt_buf[0][1] <= t[1]))
				print(srt_buf.shift().join("\t"));
			// check if multi-allelic; if not, skip the rest of part
			if (t[4].indexOf(',') < 0) {
				srt_buf.push(t); if (srt_buf.length > 1) srt_buf.sort(function(a,b) {return a[1]-b[1]});
				continue;
			}
			// parse VCF line
			var a = b8_parse_vcf2(t);
			a.sort(function(x,y){return x[0]-y[0]}); // sort by start position
			// initialize overlapping square matrix
			var co = [];
			var n_alt = t[4].split(",").length;
			for (var i = 0; i < n_alt; ++i) {
				co[i] = [];
				for (var j = 0; j < n_alt; ++j)
					co[i][j] = i == j? true : false;
			}
			// parse genotypes
			var gt = [], cnt_a = [];
			for (var i = 0; i < n_alt; ++i) cnt_a[i] = 0;
			for (var i = 9; i < t.length; ++i) {
				var m;
				if ((m = /^(\d+|\.)([\|\/])(\d+|\.)/.exec(t[i])) != null) {
					var h1 = m[1] == "."? -1 : parseInt(m[1]);
					var h2 = m[3] == "."? -1 : parseInt(m[3]);
					gt.push([h1, h2, m[2]]);
					if (h1 > 0) ++cnt_a[h1-1];
					if (h2 > 0) ++cnt_a[h2-1];
					if (h1 > 0 && h2 > 0 && h1 != h2)
						co[h1-1][h2-1] = co[h2-1][h1-1] = true;
				} else gt.push(-1, -1, '/');
			}
			for (var i = 0; i < a.length; ++i) {
				if (cnt_a[i] == 0) continue;
				// identify overlapping alleles
				var ovlp = [];
				for (var j = 0; j < a.length; ++j)
					if (i != j && cnt_a[j] > 0 && a[i][0] < a[j][1] && a[j][0] < a[i][1] && co[a[i][5]][a[j][5]])
						ovlp.push(j);
				// flag overlapping alleles
				var ovlp_flag = [];
				for (var j = 0; j < n_alt; ++j) ovlp_flag[j] = false;
				for (var j = 0; j < ovlp.length; ++j) ovlp_flag[ovlp[j]] = true;
				// update fixed fields
				t[1] = a[i][0] + 1; t[3] = a[i][3]; t[4] = a[i][4];
				if (ovlp.length > 0) t[4] += ",<M>";
				t[8] = 'GT';
				// update genotypes
				var an = a[i][5] + 1;
				for (var j = 0; j < gt.length; ++j) {
					var g = [gt[j][0], gt[j][1], gt[j][2]];
					for (var k = 0; k < 2; ++k) {
						if (gt[j][k] <= 0) continue;
						if (gt[j][k] == an) g[k] = 1;
						else if (ovlp_flag[gt[j][k] - 1]) g[k] = 2;
						else g[k] = 0;
					}
					t[j+9] = (g[0] < 0? '.' : g[0]) + g[2] + (g[1] < 0? '.' : g[1]);
				}
				srt_buf.push(t.slice()); if (srt_buf.length > 1) srt_buf.sort(function(a,b) {return a[1]-b[1]});
			}
		}
	}
	while (srt_buf.length) print(srt_buf.shift().join("\t"));

	buf.destroy();
	file.close();
}

b8_vcf2ucf(arguments);
