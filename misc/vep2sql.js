var getopt = function(args, ostr) {
	var oli; // option letter list index
	if (typeof(getopt.place) == 'undefined')
		getopt.ind = 0, getopt.arg = null, getopt.place = -1;
	if (getopt.place == -1) { // update scanning pointer
		if (getopt.ind >= args.length || args[getopt.ind].charAt(getopt.place = 0) != '-') {
			getopt.place = -1;
			return null;
		}
		if (getopt.place + 1 < args[getopt.ind].length && args[getopt.ind].charAt(++getopt.place) == '-') { // found "--"
			++getopt.ind;
			getopt.place = -1;
			return null;
		}
	}
	var optopt = args[getopt.ind].charAt(getopt.place++); // character checked for validity
	if (optopt == ':' || (oli = ostr.indexOf(optopt)) < 0) {
		if (optopt == '-') return null; //  if the user didn't specify '-' as an option, assume it means null.
		if (getopt.place < 0) ++getopt.ind;
		return '?';
	}
	if (oli+1 >= ostr.length || ostr.charAt(++oli) != ':') { // don't need argument
		getopt.arg = null;
		if (getopt.place < 0 || getopt.place >= args[getopt.ind].length) ++getopt.ind, getopt.place = -1;
	} else { // need an argument
		if (getopt.place >= 0 && getopt.place < args[getopt.ind].length)
			getopt.arg = args[getopt.ind].substr(getopt.place);
		else if (args.length <= ++getopt.ind) { // no arg
			getopt.place = -1;
			if (ostr.length > 0 && ostr.charAt(0) == ':') return ':';
			return '?';
		} else getopt.arg = args[getopt.ind]; // white space
		getopt.place = -1;
		++getopt.ind;
	}
	return optopt;
}

// recommended VEP command line:
//   ./variant_effect_predictor.pl -i input.vcf -o output.txt --chr 1 --offline --pick --cache --fork 16 --everything --force_overwrite --quiet

var schema = "\
CREATE TABLE Variant (\n\
  vid         TEXT,    -- unique variant ID in the format of chr:pos:rlen:base\n\
  chrom       TEXT,    -- contig name\n\
  bin         INTEGER, -- BAM bin number\n\
  chromStart  INTEGER, -- contig start, 0-based (BED-like)\n\
  chromEnd    INTEGER, -- contig end, BED-like\n\
  bases       TEXT,    -- allele sequence\n\
  impact      INTEGER, -- 0=MODIFER, 1=LOW, 2=MODERATE, 3=HIGH (most significant only)\n\
  effect      TEXT,    -- Sequence Ontology term (most significant only)\n\
  gene        TEXT,    -- stable identifer of gene symbol\n\
  biotype     TEXT,    -- biotype of transcript or regulatory feature\n\
  cdsPos      INTEGER, -- position on CDS\n\
  featID      TEXT,    -- typically the transcript ID\n\
  distance    INTEGER, -- shortest distance from variant to transcript\n\
  strand      INTEGER, -- strand; 1 or -1\n\
  codonChg    TEXT,    -- codon change\n\
  ccds        TEXT,    -- CCDS ID\n\
  sift        TEXT,    -- SIFT effect\n\
  polyphen    TEXT,    -- PolyPhen effect\n\
  PRIMARY KEY (vid)\n\
);\n\
";

var severity = [
	"transcript_ablation",
	"splice_acceptor_variant",
	"splice_donor_variant",
	"stop_gained",
	"frameshift_variant",
	"stop_lost",
	"start_lost",
	"transcript_amplification",
	"inframe_insertion",
	"inframe_deletion",
	"missense_variant",
	"protein_altering_variant",
	"splice_region_variant",
	"incomplete_terminal_codon_variant",
	"stop_retained_variant",
	"synonymous_variant",
	"coding_sequence_variant",
	"mature_miRNA_variant",
	"5_prime_UTR_variant",
	"3_prime_UTR_variant",
	"non_coding_transcript_exon_variant",
	"intron_variant",
	"NMD_transcript_variant",
	"non_coding_transcript_variant",
	"upstream_gene_variant",
	"downstream_gene_variant",
	"TFBS_ablation",
	"TFBS_amplification",
	"TF_binding_site_variant",
	"regulatory_region_ablation",
	"regulatory_region_amplification",
	"feature_elongation",
	"regulatory_region_variant",
	"feature_truncation",
	"intergenic_variant"
];

var impact_str2int = {"MODIFIER":0, "LOW":1, "MODERATE":2, "HIGH":3};

var severity_rank = {};
for (var i = 0; i < severity.length; ++i)
	severity_rank[severity[i]] = i;

function reg2bin(start, end)
{
    --end;
    if (start>>14 == end>>14) return ((1<<15)-1)/7 + (start>>14);
    if (start>>17 == end>>17) return ((1<<12)-1)/7 + (start>>17);
    if (start>>20 == end>>20) return ((1<<9)-1)/7  + (start>>20);
    if (start>>23 == end>>23) return ((1<<6)-1)/7  + (start>>23);
    if (start>>26 == end>>26) return ((1<<3)-1)/7  + (start>>26);
    return 0;
}

var c, is_quiet = false, is_fmf = false, create_tbl = false;
while ((c = getopt(arguments, "hqfc")) != null) {
	if (c == 'q') is_quiet = true;
	else if (c == 'f') is_fmf = true;
	else if (c == 'c') create_tbl = true;
	else if (c == 'h') {
		print("Usage: k8 vep2sql.js [-fch] <vep-out.txt>");
		print("Options:");
		print("  -q      quiet mode");
		print("  -f      output in the FMF format (default is SQL)");
		print("  -c      output SQL table schema");
		print("");
		print("Convert VEP native output to SQLite schema or the FMF format. It has been");
		print("tested on VEP v80 on command line options:");
		print("  -i input.vcf -o vep-out.txt --offline --pick --cache --everything");
		exit(1);
	}
}

var file = arguments.length > getopt.ind? new File(arguments[getopt.ind]) : new File();
var buf = new Bytes();
var lineno = 0;

if (!is_fmf) {
	if (create_tbl) print(schema);
	print("BEGIN TRANSACTION;");
}
while (file.readline(buf) >= 0) {
	++lineno;
	if (buf[0] == 35) continue; // a # line
	var m, t = buf.toString().split("\t");
	var info = [];

	// generate vid
	var chr, start, end;
	if ((m = /^(\S+):(\d+)(-(\d+))?/.exec(t[1])) != null) {
		chr = m[1];
		start = parseInt(m[2]);
		end = m[3]? parseInt(m[4]) : start;
	} else {
		if (!is_quiet)
			warn("["+lineno+"] failed to parse contig name and position");
		continue;
	}
	var bases, rlen, is_sym = false;
	if (t[2] == '-') { // deletion
		bases = "", rlen = end - start + 1;
	} else if (start == end && t[2].length == 1) { // SNP
		bases = t[2], rlen = 1;
	} else if (end - start == 1) { // insertion
		if (t[2] == 'insertion') {
			bases = "<INS>";
			is_sym = true;
		} else bases = t[2];
		rlen = 0, start = end;
	} else {
		bases = "<" + t[2] + ">";
		rlen = end - start + 1;
		is_sym = true;
	}
	var vid = [chr, start, rlen, bases].join(":");
	if (is_sym) {
		if (!is_quiet)
			warn("["+lineno+"] skipped symbolic allele " + vid);
		continue;
	}

	// write positions
	if (!is_fmf) {
		info.push(["vid", "Z", vid]);
		info.push(["chrom", "Z", chr]);
		info.push(["bin", "i", reg2bin(start - 1, start - 1 + rlen)]);
		info.push(["chromStart", "i", start - 1]);
		info.push(["chromEnd", "i", start - 1 + rlen]);
		info.push(["bases", "Z", bases]);
	}

	// parse other fields
	if ((m = /SYMBOL=([^\s;]+);SYMBOL_SOURCE=HGNC;.*BIOTYPE=([^\s;]+)/.exec(t[13])) != null) {
		info.push(["gene", "Z", m[1]]);
		info.push(["biotype", "Z", m[2]]);
	}
	if ((m = /IMPACT=([^\s;]+)/.exec(t[13])) != null) {
		if (impact_str2int[m[1]] != null)
			info.push(["impact", "i", impact_str2int[m[1]]]);
	}
	if (t[6].indexOf(",") > 0) {
		var u = t[6].split(","), max = -1, max_eff = null;
		for (var j = 0; j < u.length; ++j) {
			var r = severity_rank[u[j]];
			if (r == null) {
				if (!is_quiet)
					warn("["+lineno+"] unknown effect " + u[j]);
				continue;
			}
			if (r > max) max = r, max_eff = u[j];
		}
		info.push(["effect", "Z", max_eff]);
	} else info.push(["effect", "Z", t[6]]);
	if (t[4] != '-' && t[5] != '-')
		info.push(["featID", "Z", t[4]]);
	if (t[8] != '-') info.push(["cdsPos", "i", t[8]]);
	if ((m = /DISTANCE=(\d+);STRAND=(-?\d+)/.exec(t[13])) != null) {
		info.push(["distance", "i", m[1]]);
		info.push(["strand", "i", m[2]]);
	}
	if ((m = /;CCDS=([^\s;]+)/.exec(t[13])) != null) info.push(["ccds", "Z", m[1]]);
	if (t[10] != '-') info.push(["codonChg", "Z", t[11]]);
	if ((m = /;SIFT=([^\s;()]+)\(([\d.]+)\)/.exec(t[13])) != null) info.push(["sift", "Z", m[1]]);
	if ((m = /;PolyPhen=([^\s;()]+)\(([\d.]+)\)/.exec(t[13])) != null) info.push(["polyphen", "Z", m[1]]);

	// generate SQL or FMF line
	if (!is_fmf) {
		var key = [], val = [];
		for (var i = 0; i < info.length; ++i) {
			key.push(info[i][0]);
			if (info[i][1] != 'Z') val.push(info[i][2]);
			else val.push("'" + info[i][2] + "'");
		}
		print("INSERT INTO Variant (" + key.join(",") + ") VALUES (" + val.join(",") + ");");
	} else {
		var out = [vid];
		for (var i = 0; i < info.length; ++i)
			out.push(info[i].join(":"));
		print(out.join("\t"));
	}
}
if (!is_fmf) {
	print("END TRANSACTION;\n");
	print("CREATE INDEX idx_gene  ON Variant (gene);");
	print("CREATE INDEX idx_bin   ON Variant (chrom, bin);");
	print("CREATE INDEX idx_start ON Variant (chrom, chromStart);");
}

buf.destroy();
file.close();
