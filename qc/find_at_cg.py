import subprocess

path="/scratch/cancgene/plaw/CLL_spanish/plink"
with open("/scratch/cancgene/plaw/CLL/snp genotypes/replication_new") as repl_file:
    for line in repl_file:
        rsid, chrom, pos = line.rstrip().split()
        print(rsid)
        fname="{}/caco/{}_chr{}_caco-merge.missnp".format(path,rsid, chrom)
        fout=fname+".clean"
        fout2=fname+".exclude"
        with open(fname) as missnpfile, open(fout, "w") as ofile, open(fout2, "w") as ofile2:
            for missnp in missnpfile:
                missnp=missnp.rstrip()
                #print("#",missnp)
                res=subprocess.check_output("grep -m 1 -w {} {}/separate/{}_chr{}_controls.bim".format(missnp, path, rsid, chrom), shell=True).decode("utf-8")
                tokens=res.split()
                a1=tokens[-2]
                a2=tokens[-1]
                if a1=="0" or a2=="0":
                    ofile2.write(missnp+"\n")
                elif (a1=="C" and a2=="G") or (a1=="G" and a2=="C") or (a1=="A" and a2=="T") or (a1=="T" and a2=="A"):
                    ofile2.write(missnp+"\n")
                else:
                    ofile.write(missnp+"\n")

