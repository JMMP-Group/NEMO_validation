
#!/bin/bash
echo "Bash version ${BASH_VERSION}..."
MOD=$1
for (( month=1; month<13; month++ ))
do 
     sbatch  merge.sbatch  ${MOD} ${month}
done
