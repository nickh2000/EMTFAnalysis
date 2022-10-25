

for i in $(seq 1 $2); do
    if [ $i -eq $2 ];
    then
        echo Last!

        python $1 -n $2 -i $(( $i - 1 )) &
    else
        echo Not Last!
        python $1 -n $2 -i $(( $i - 1 )) &
    fi
    echo Running job: $i of $2
    pids[${i}]=$!
done