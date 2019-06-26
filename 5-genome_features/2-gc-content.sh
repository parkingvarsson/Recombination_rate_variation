#! /bin/bash -l
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00

##### Script for calculating gc-content across the genome from a fasta file #####
# $1 should be the reference fasta -file
##### Edited from: https://www.biostars.org/p/70167/#70172 #####

awk ' \
BEGIN { \
    FS=""; \
    cg=0; \
    t=0; \
    gap=0; \
    cw=0; \
    pw=0; \
    step = 50000; \
} \
{ \
    if ($1 != ">") { \
        for (i = 1; i <= NF; i++) { \
            if ($i ~ /[ACTGNactgn]/) { \
                t++;
            } \
            if ($i ~ /[Nn]/) { \
                gap++;
            } \
            if ($i ~ /[CGcg]/) { \
                cg++;
            } \
            if (t == step) { \
                if (gap/t > 0.8) { \
                    pw = cw; \
                    cw = cw+t; \
                    print h"\t"pw"\t"cw"\t""NA""\t""NA""\t""NA""\t""NA";
                    cg = 0; \
                    t = 0; \
                    gap = 0; \
                } \
                else { \
                    pw = cw; \
                    cw = cw+t; \
                    print h"\t"pw"\t"cw"\t"cg"\t"t-gap"\t"gap"\t"(cg/(t-gap)); \
                    cg = 0; \
                    t = 0; \
                    gap = 0; \
                } \
            } \
        } \
    } \
    else { \
        if (t > 0) { \
            if (gap/t > 0.8) { \
                pw = cw; \
                cw = cw+t; \
                print h"\t"pw"\t"cw"\t""NA""\t""NA""\t""NA""\t""NA";
                cg = 0; \
                t = 0; \
                gap = 0; \
                pw=0; \
                cw=0; \
            } \
            else { \
                pw = cw; \
                cw = cw+t; \
                print h"\t"pw"\t"cw"\t"cg"\t"t-gap"\t"gap"\t"(cg/(t-gap)); \
                cg = 0; \
                t = 0; \
                gap = 0; \
                pw=0; \
                cw=0; \
            } \
        } \
        h = substr($0,2); \
    } \
} \
END { \
    if (t > 0) { \
        pw = cw; \
        cw = cw+t; \
        if (gap/t > 0.8) { \
            print h"\t"pw"\t"cw"\t""NA""\t""NA""\t""NA""\t""NA";
        } \
        else { \
           print h"\t"pw"\t"cw"\t"cg"\t"t-gap"\t"gap"\t"(cg/(t-gap)); \
           cg = 0; \
           t = 0; \
        } \
    } \
}' $1 > gc-content.txt
