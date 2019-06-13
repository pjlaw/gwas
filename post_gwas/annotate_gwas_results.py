# -*- coding: utf-8 -*-
"""
Script to annotate CRC GWAS results
Updated from Sara's version
Feb 2019

@author: plaw
"""
import sys
import sqlite3
from math import exp
conn = sqlite3.connect(":memory:")
c=conn.cursor()

c.execute('''CREATE TABLE cytoband_dat
                 (chromosome text, start integer, end integer, cytoband text)''')

#user output options
fname=sys.argv[1] #output file name
data_dir=sys.argv[2]
outfile="{}/{}_annotated_results.txt".format(data_dir, fname)

ld_threshold=0.1; # annotates snps as in LD with known variants if r2 is greater than this value

clump_file="{}/plink_res/{}_clump_all.clumped".format(data_dir, fname)
meta_results="{}/{}.txt".format(data_dir, fname)
header_file="{}/head.txt".format(data_dir)

#known data constants
known_snps_file="/scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/US_meta/known_snps_eur_plink.txt"
#known_snps_file="/scratch/DGE/MOPOPGEN/plaw/CRC_GWAS/scotland_meta/output/crc_published_snps.txt"
ld_data="{}/plink_res/LD_all.ld".format(data_dir)
cytoband_file="/scratch/DGE/MOPOPGEN/plaw/reference_data/cytoBand-hg19.txt"

known_snps={}
chromosome_data={}
print("parsing known snps file")
try:
    with open(known_snps_file) as ifile:
        for line in ifile:
            rsid, chrom, pos = line.rstrip().split()
            known_snps[rsid]=[chrom, pos]
            if chrom in chromosome_data.keys():
                chromosome_data[chrom].append(int(pos))
            else:
                chromosome_data[chrom]=[int(pos)]
except FileNotFoundError:
    sys.exit("known snps file not found")

print("parsing cytoband file")
try:
    with open(cytoband_file) as ifile:
        for line in ifile:
            bits = line.rstrip().split()
            c.execute("INSERT INTO cytoband_dat VALUES (?, ?, ?, ?)", [bits[0][3:], int(bits[1])+1, int(bits[2]), bits[3]])
except FileNotFoundError:
    sys.exit("cytoband file not found")

c.execute("create index cyto_idx on cytoband_dat(chromosome, start, end)")
conn.commit()

print("parsing known snp LD file")
gwas_ld={}
try:
    with open(ld_data) as ifile:
        for line in ifile:
            bits = line.strip().split()
            r2=float(bits[6])
            dprime=bits[7]
            snp_b=bits[5]
            snp_a=bits[2]
            if r2>ld_threshold:
                if snp_b not in known_snps.keys():
                    #ignore if in ld with known snp
                    if snp_b in gwas_ld.keys():
                        gwas_ld[snp_b]=gwas_ld[snp_b]+",{}({},{})".format(snp_a, r2, dprime)
                    else:
                        gwas_ld[snp_b]="{}({},{})".format(snp_a, r2, dprime)

except FileNotFoundError:
    sys.exit("LD data file not found")


count_clump=1
snp_hash={}
lead_snp_hash={}
number_supports={}
print("parsing clumped data file")
try:
    with open(clump_file) as ifile:
        for line in ifile:
            bits = line.strip().split()
            if len(bits)==0:
                #skip out the blank lines from the cat
                continue
            lead_snp=bits[2]
            other_snps=bits[11].split("(1),")
            count_str=str(count_clump)
            lead_snp_hash[lead_snp]=count_str #rsid=clump_id
            snp_hash[lead_snp]=[count_str]

            for snppie in other_snps:
                if snppie=="":
                    continue
                if snppie in snp_hash.keys():
                    snp_hash[snppie].append(count_str)
                else:
                    snp_hash[snppie]=[count_str]

            number_supports[count_str]=str(len(other_snps))
            count_clump+=1
except FileNotFoundError:
    sys.exit("clump file not found")


header=""
try:
    with open(header_file) as ifile:
        header=ifile.readline()
except FileNotFoundError:
    sys.exit("header file not found")

header_dict=dict((value, key) for (key, value) in enumerate(header.strip().split()))

print("generating output")
try:
    with open(meta_results) as ifile, open(outfile, "w") as ofile:
        ofile.write("rsid cytoband lead known distance_to_known clump support ld_to_known studies RA/other RAF OR(95%CI) {}".format(header))

        for line in ifile:
            data=line.rstrip().split()
            rsid=data[1]
            chrom=data[0]
            pos=data[2]
            outline=[rsid]

            #add the cytoband
            try:
                cyto_res = c.execute("SELECT cytoband FROM cytoband_dat where chromosome={} and start<{} and end>{}".format(chrom, pos, pos)).fetchone()[0]
                outline.append("{}{}".format(chrom, cyto_res))
            except TypeError:
                outline.append("-")


            #check if this is the lead snp, ie top p-value
            if rsid in lead_snp_hash.keys():
                outline.append("LEAD")
            else:
                outline.append("-")

            #if this a known snp
            if rsid in known_snps.keys():
                outline.append("KNOWN")
            else:
                outline.append("-")

            #calculate the distance to the closest known snp
            distance=sys.maxsize
            try:
                for known_position in chromosome_data[chrom]:
                    diff=abs(known_position-int(pos))
                    if diff<distance:
                        distance=diff
                if distance==sys.maxsize:
                    outline.append("-")
                else:
                    outline.append(str(distance))
            except KeyError:
                outline.append("-")

            #output clump and support data
            if rsid in snp_hash.keys():
                this_ld_clumps=snp_hash[rsid]
                outline.append("clump{}".format(";".join(this_ld_clumps)))
                this_supports=[]
                for ld_clump in this_ld_clumps:
                    this_supports.append(number_supports[ld_clump])
                outline.append(";".join(this_supports))
            else:
                outline+=["-","-"]

            #if this snp in ld with a known snp
            try:
                outline.append(gwas_ld[rsid])
            except KeyError:
                outline.append("-")

            #count the number of studies this snp is in
            study_count=0
            startpos=header_dict["cohort_1_p"]
            for i in data[startpos::5]:#first p-value at index X, every 5th column
                if float(i)>0:
                    study_count+=1
            outline.append(str(study_count))

            #calculate risk alleles
            beta=float(data[header_dict["beta"]])
            se=float(data[header_dict["se"]])
            ra_str=""
            raf=0
            if beta>0:
                #allele_B is the risk allele
                ra_str="{}/{}".format(data[4], data[3])
                raf=data[6]
            else:
                ra_str="{}/{}".format(data[3], data[4])
                raf=str(1-float(data[6]))
            outline+=[ra_str, raf]

            # calulcate odds ratio
            beta=abs(beta)
            or_str="{:.2f}({:.2f}-{:.2f})".format(exp(beta), exp(beta-1.96*se), exp(beta+1.96*se))
            outline.append(or_str)

            #add the meta data
            outline=outline+data

            ofile.write(" ".join(outline)+"\n")
except FileNotFoundError:
    sys.exit("result file not found")


conn.close()
print("done")
