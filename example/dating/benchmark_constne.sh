#for SEED in {1..20}; do
for SEED in {21..30}; do
  echo "bash simulate_relate_constne.sh $SEED &>log.$SEED"
done >jobfile.constne
parallel -j 12 <jobfile.constne
