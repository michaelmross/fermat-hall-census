\\ skeleton_ranks.gp -- unconditional Mordell-Weil rank census of the curve
\\ classes arising from the {2,3,m} anchor ensemble.
\\
\\ Input:  skeleton_classes.txt, one sixth-power-free integer k0 per line,
\\         written by `python coprimality_pipeline.py skeleton`.
\\ Output: skeleton_ranks.txt, one line per class:
\\
\\             k0|rank_lower|rank_upper|torsion_order
\\
\\         The rank is determined unconditionally when the bounds coincide;
\\         ERR marks a class PARI could not resolve.
\\
\\ Usage:  gp -q skeleton_ranks.gp
\\ Tested with PARI/GP 2.15.4 (ellrank was introduced in 2.14).
\\
\\ NOTE: gp parses script files one line at a time, so every construct below
\\ is kept on a single line. Do not reflow.

default(parisize, "512M");
infile = "skeleton_classes.txt";
outfile = "skeleton_ranks.txt";
classes = readvec(infile);
n = #classes;
print("classes to process: ", n);
write(outfile, "");
doclass(k) = iferr(E = ellinit([0, k]); r = ellrank(E); t = elltors(E)[1]; write1(outfile, Str(k, "|", r[1], "|", r[2], "|", t, "\n")); if(r[1] == r[2], 1, 0), err, write1(outfile, Str(k, "|ERR|ERR|ERR\n")); -1);
det = 0; amb = 0; bad = 0;
for(i = 1, n, s = doclass(classes[i]); if(s == 1, det++, if(s == 0, amb++, bad++)); if(i % 200 == 0, print("  ", i, "/", n)));
print("rank determined (bounds coincide): ", det);
print("rank bounds ambiguous:             ", amb);
print("errors:                            ", bad);
print("wrote ", outfile);
quit
