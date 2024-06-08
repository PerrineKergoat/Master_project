import pyslim
import tskit
import msprime
import sys

job_id = sys.argv[1]
rep_nb = sys.argv[2]
recomb_rate = sys.argv[3]
gen_nb = sys.argv[4]
nb_pop = sys.argv[5]
path_data = sys.argv[6]

### Recapitation
orig_ts = tskit.load(f"{path_data}Simu_{gen_nb}_{job_id}_{rep_nb}.trees")

print("File opened")
print(f"Replicates : {rep_nb}")

rts = pyslim.recapitate(orig_ts,
            recombination_rate=recomb_rate,
            ancestral_Ne=1000*int(nb_pop)).simplify()

rts = pyslim.generate_nucleotides(rts)
rts = pyslim.convert_alleles(rts)

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())

print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

### Add neutral mutations
for var in rts.variants():
    print(var.site.position, var.alleles, var.genotypes, sep="\t")

rts = msprime.sim_mutations(rts, rate=1e-7, model = msprime.BinaryMutationModel())

### Create VCF output
with open(f"{path_data}Recap_{gen_nb}_{job_id}_{rep_nb}.vcf", "w") as simuvcf:
    rts.write_vcf(simuvcf, position_transform = "legacy")

print("Done")
