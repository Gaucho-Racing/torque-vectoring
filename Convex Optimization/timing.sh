#!/bin/bash

# Configuration
program="./NewtonReworked"   # Change this to your executable
runs=20                    # Number of times to run
total=0
max=0.
min=10000.
iter=0
yawMoment=0
totalTorque=0
steeringAngle=0
bad=0
echo "Running $program $runs times..."
for((k = -400; k <= 800; k += 20)) do
	for ((j = -30; j <= 30; j += 1)); do
		for ((i = -500 ; i <= 500; i+=20)); do
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
				yawMoment=$i
				steeringAngle=$j
				totalTorque=$k
			fi

		    else
			echo "Warning: Non-numeric output '$output' on run $j $i $k"
			bad=$((bad+1))
		    fi
		done
#		echo "Angle $j Done"
	done
	echo "Total Torque $k Done"
done

# Calculate average
average=$(echo "scale=4; $total / $iter" | bc)
badPer=$(echo "scale=4; $bad / $iter" | bc)

echo "Average output: $average"
echo "Bad Percentage: $badPer"
echo "Max Time: $max"
echo "Min Time: $min"
echo "$steeringAngle $yawMoment $totalTorque"
