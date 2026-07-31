default(parisize, "512M");
v = readvec("skeleton_classes.txt");
n = #v; nsha = 0; shalist = List();
for(i = 1, n, E = ellinit([0, v[i]]); r = ellrank(E); if(r[3] > 0, nsha++; listput(shalist, [v[i], r[3], r[2]])));
print("classes with nontrivial Sha[2] indicator: ", nsha, " of ", n);
for(i = 1, #shalist, print("  k0=", shalist[i][1], "  s=", shalist[i][2], "  rank<=", shalist[i][3]));
quit
