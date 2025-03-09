#!/bin/bash
# Generate a random using `shuf` command
echo "How many random numbers do you want to generate?:"
read number

echo "nmax=?"
read nmax

echo
echo
#Print the generated random numbers
entries=($(shuf -i 1-${nmax} -n ${number}))

rm -rf rand_num.dat
for entry in "${entries[@]}"; do
    echo "$entry" >> rand_num.dat
done

