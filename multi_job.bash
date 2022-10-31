#This program exists to parallelize other python scripts by splitting up input files to scripts into chunks
#Argument #1: program to paralellize
#Argument #2: number of batches
#Iterate through number of batches given (second argument)
for i in $(seq 1 $2); do
    #Note the last job submitted
    if [ $i -eq $2 ];
    then
        echo Last!
        #Python scripts are programmed to take argument -n (number of jobs) and -i (index of batch), which are both needed to grab a chunk of the input files
        python $1 -n $2 -i $(( $i - 1 )) &
    else
        echo Not Last!
        python $1 -n $2 -i $(( $i - 1 )) &
    fi
    echo Running job: $i of $2

    
    pids[${i}]=$!
done

#This was used before to mark which job-processes were still running and recombine them at the end with recombine.py, since depreciated
#Can be pretty easily re-implemented

# for pid in ${pids[*]}; do
#     wait $pid
# done