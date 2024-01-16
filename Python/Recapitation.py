import pyslim 
import tskit 
import msprime
import sys

replicateN = sys.argv[1]

### Recapitation 

orig_ts = tskit.load(f"/work/FAC/FBM/DEE/jgoudet/default/pkergoat/2pop_1sel_500snps_100overlap/data/TestFST_{replicateN}.trees")

print("File opened")

rts = pyslim.recapitate(orig_ts,
            recombination_rate=1e-8,
            ancestral_Ne=2000).simplify()

rts = pyslim.generate_nucleotides(rts)
rts = pyslim.convert_alleles(rts)

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())

print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

### Add mutations 

for var in rts.variants():
    print(var.site.position, var.alleles, var.genotypes, sep="\t")

#model = msprime.SLiMMutationModel(type=1)
rts = msprime.sim_mutations(rts, rate=1e-7) #, random_seed=5678) #, model=model)
#print("test")

#for var in rts.variants():
    #print(var.site.position, var.alleles, var.genotypes, sep="\t")

 

### Create VCF output

with open(f"/work/FAC/FBM/DEE/jgoudet/default/pkergoat/2pop_1sel_500snps_100overlap/data/Simu_testFST_{replicateN}.vcf", "w") as simuvcf:
    rts.write_vcf(simuvcf, position_transform = "legacy")

print("done")
