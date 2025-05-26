#!/bin/bash

# Configuration
program="./NewtonReworked 0 300 430"   # Change this to your executable
runs=20                    # Number of times to run
total=0
max=0.
min=10000.
iter=0
echo "Running $program $runs times..."
for((k = -400; k <= 800; k += 100)) do
	for ((j = -30; j <= 30; j++)); do
		for ((i = -500 ; i <= 500; i+=100)); do
		    program="./NewtonReworked $j $i $k"
		    output=$($program)
	#	    echo "Run $i: $output"

		    # Make sure output is numeric
		    if [[ $output =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
			total=$(echo "$total + $output" | bc)
			iter=$((iter + 1))	
			if [[ $output < $min ]]; then
				min=$output
			fi
			if [[ $output > $max ]]; then
				max=$output
			fi

		    else
			echo "Warning: Non-numeric output '$output' on run $j $i $k"
		    fi
		done
#		echo "Angle $j Done"
	done
	echo "Total Torque $k Done"
done

# Calculate average
average=$(echo "scale=4; $total / $iter" | bc)

echo "Average output: $average"
echo "Max Time: $max"
echo "Min Time: $min"
