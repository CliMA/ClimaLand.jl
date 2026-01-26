#!/bin/bash

julia benchmark/send_comment_to_pr.jl -o $org -r $repo -p $pullrequest -c '**Starting benchmarks!**'  

if julia benchmark/run_benchmarks.jl $1 ; then
    julia benchmark/send_comment_to_pr.jl -o $org -r $repo -p $pullrequest -g
else
    julia benchmark/send_comment_to_pr.jl -o $org -r $repo -p $pullrequest -c "**An error has occured while running the benchmarks in file $1** "
fi
