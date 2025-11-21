nfold="fold-restart"

mkdir $nfold
cp step-1/input.yaml $nfold
cp step-1/job_fit.sh $nfold

sed -i 's|pacemaker input.yaml|pacemaker input.yaml -p ../step-1/output_potential.yaml|g' $nfold/job_fit.sh
sed -i 's/kappa: 0.95/kappa: 0.30/g' $nfold/input.yaml
sed -i 's|filename: |filename: ../step-1/|g' $nfold/input.yaml
