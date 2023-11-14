set -e

#binaries from website
RELATE_DIR=`pwd`/../../../../repo/relate
#git clone https://github.com/leospeidel/relate_lib.git
#cd relate_lib && mkdir -p build && cd build
#cmake .. && make && cd ../..
RELATE_LIB=`pwd`/../../../../relate_lib 

SEED="2000"
ID="chr1"
SEQLEN=10000000
SAMPLES=1000 #diploid
MCMC_DRAWS=1000
QUAD_ORDER=10000 #importance samples
NUM_THREADS=12

WORK_DIR=`pwd`/sim_varne2/$SEED
IN_DIR=$WORK_DIR/relate_inputs
OUT_DIR=$WORK_DIR/relate_outputs
BENCH_DIR=$WORK_DIR/benchmarks
mkdir -p $WORK_DIR
mkdir -p $IN_DIR
mkdir -p $OUT_DIR
mkdir -p $BENCH_DIR


# --- SIMULATE DATA --- #

cd $WORK_DIR
#echo "SEED:$SEED ID:$ID SEQLEN:$SEQLEN SAMPLES:$SAMPLES" >$ID.info
#printf '
#import stdpopsim
#import numpy as np
#from sys import stdout
#seed = %s
#id = "%s"
#seqlen = %s
#samples = %s
#homsap = stdpopsim.get_species("HomSap")
#zigzag = homsap.get_demographic_model("Zigzag_1S14")
#engine = stdpopsim.get_engine("msprime")
#genmap = homsap.get_contig("chr1", genetic_map="HapMapII_GRCh38").recombination_map
#endpnt = genmap.position[1]
#chrom1 = homsap.get_contig(
#  "chr1", genetic_map="HapMapII_GRCh38", 
#  left=endpnt, right=endpnt+seqlen,
#)
#muout = open(f"{id}.mu", "w")
#print(chrom1.mutation_rate, file=muout)
#muout.close()
#ts = engine.simulate(
#  contig=chrom1, demographic_model=zigzag, 
#  samples={"generic":samples}, seed=seed,
#)
#ts.dump(f"{id}.trees")
#ts.write_vcf(
#  open(f"{id}.vcf", "w"), 
#  contig_id=id,
#  individual_names=[f"ZIGZAG{i}" for i in range(samples)],
#  position_transform=lambda p: [x-1 for x in p],
#)
#neout = open(f"{id}.ne", "w")
#print(int(0.5 * ts.diversity() / chrom1.mutation_rate), file=neout)
#neout.close()
#
#lnout = open(f"{id}.seqlen", "w")
#print(int(ts.sequence_length), file=lnout)
#lnout.close()
#
#dbg = zigzag.model.debug()
#steps = np.linspace(0, 40000, 20)
#coalrate, _ = dbg.coalescence_rate_trajectory(lineages={"generic":2}, steps=steps)
#coalout = open(f"{id}.coal", "w")
#print("0", file=coalout)
#print(" ".join([str(int(x)) for x in steps]), file=coalout)
#print("0 0", " ".join([str(x) for x in coalrate]), file=coalout)
#coalout.close()
#recout = open(f"{id}.hapmap", "w")
#genmap = chrom1.recombination_map
#breaks = genmap.position
#cummas = genmap.get_cumulative_mass(breaks)
#rates = np.concatenate([genmap.rate * 1e8, [0.0]])
#print("Position(bp) Rate(cM/Mb) Map(cM)", file=recout)
#for b, c, r in zip(breaks, cummas, rates):
#  print(int(b), r, c, file=recout)
#recout.close()
#' $SEED $ID $SEQLEN $SAMPLES | python3 && gzip -f $ID.vcf

# ----------- INFER TREES W TSINFER ----------- #

#cd $WORK_DIR
#printf "
#import tsinfer
#import tskit
#
#ts_true = '%s'
#ts_infr = '%s'
#ts_true = tskit.load(ts_true)
#sample_data = tsinfer.SampleData.from_tree_sequence(ts_true)
#ts = tsinfer.infer(sample_data).simplify(filter_sites=False)
#ts.dump(ts_infr)
#" $WORK_DIR/chr1.trees $WORK_DIR/chr1.infer.trees | python

cd $WORK_DIR
MU=`cat $ID.mu`
NE=`cat $ID.ne`
SEQLEN=`cat $ID.seqlen`


# --- INFER TREES VIA RELATE --- #

cd $IN_DIR

##reassign reference alleles
#bioawk -cvcf -H '{$alt="G"; $ref="A"; print}' ../$ID.vcf.gz >$ID.vcf
#
##dump into haplotypes
#$RELATE_DIR/bin/RelateFileFormats \
#  --mode ConvertFromVcf --haps $ID.haps \
#  --sample $ID.samples -i $ID
#
##hapmap
#cp $WORK_DIR/$ID.hapmap .
#
##coal rate
#cp $WORK_DIR/$ID.coal .
#
##population labels
#awk 'NR==1 {print "sample\tpopulation\tgroup\tsex"} NR>2{pop=$1; sub(/[0-9]+/, "", pop); print $1"\t"pop"\t"pop"\t2"}' $ID.samples >$ID.poplabels
#
##ancestral states
#printf '
#cat(">%s", "\n", sep="")
#cat(rep("A", %s), "\n", sep="")
#' $ID $SEQLEN | R --slave >$ID.anc.fa
#
##accessibility mask
#printf '
#cat(">%s", "\n", sep="")
#cat(rep("P", %s), "\n", sep="")
#' $ID $SEQLEN | R --slave >$ID.mask.fa
#
##infer trees
#$RELATE_DIR/scripts/PrepareInputFiles/PrepareInputFiles.sh \
#  --haps $ID.haps \
#  --sample $ID.samples \
#  --poplabels $ID.poplabels \
#  --ancestor $ID.anc.fa \
#  --mask $ID.mask.fa \
#  -o $ID

## this failed with the beta < 0 error
#$RELATE_DIR/bin/Relate --mode All -m $MU -N $NE \
#  --haps $ID.haps.gz --sample $ID.sample.gz --annot $ID.annot --seed $SEED \
#  --coal $ID.coal \
#  --variational --order $QUAD_ORDER \
#  --dist $ID.dist.gz --map $ID.hapmap -o $ID

#$RELATE_DIR/bin/Relate --mode All -m $MU -N $NE \
#  --haps $ID.haps.gz --sample $ID.sample.gz --annot $ID.annot --seed $SEED \
#  --coal $ID.coal \
#  --dist $ID.dist.gz --map $ID.hapmap -o ${ID}_mcmc

# --- ESTIMATE BRANCH LENGTHS USING MCMC --- #

## --- DUMP TRUE TOPOLOGIES
#cd $WORK_DIR
#$RELATE_LIB/bin/Convert --mode ConvertFromTreeSequence \
#  -i $WORK_DIR/$ID.trees --anc $OUT_DIR/${ID}_true.anc --mut $OUT_DIR/${ID}_true.mut
#
## --- ESTIMATE BRANCH LENGTHS USING MCMC ON TRUE TOPOLOGIES --- #
#cd $WORK_DIR
#$RELATE_DIR/scripts/SampleBranchLengths/ReEstimateBranchLengths.sh \
#  -i $OUT_DIR/${ID}_true \
#  -o $OUT_DIR/${ID}_true.sample \
#  -m $MU \
#  --coal $WORK_DIR/$ID.coal \
#  --seed $SEED --threads 10

# --- CONVERT TO TREE SEQ --- #
#cd $WORK_DIR
#$RELATE_LIB/bin/Convert --mode ConvertToTreeSequence \
#  -o $OUT_DIR/${ID}_true.sample --anc $OUT_DIR/${ID}_true.sample.anc --mut $OUT_DIR/${ID}_true.sample.mut

#cd $WORK_DIR
#$RELATE_LIB/bin/Convert --mode ConvertToTreeSequence \
#  -o $OUT_DIR/${ID}_mcmc --anc $IN_DIR/${ID}_mcmc.anc --mut $IN_DIR/${ID}_mcmc.mut
