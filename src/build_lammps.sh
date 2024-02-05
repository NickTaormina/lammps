#!/bin/bash

# Step 1: Build 'yes-USER-CAC'
echo "Building yes-USER-CAC..."
make yes-USER-CAC

# Check if the build was successful
if [ $? -ne 0 ]; then
    echo "Error: Build of 'yes-USER-CAC' failed."
    exit 1
fi

# Step 2: Build 'mpi'
echo "Building mpi..."
make mpi

# Check if the build was successful
if [ $? -ne 0 ]; then
    echo "Error: Build of 'mpi' failed."
    exit 1
fi

# Step 3: Copy 'lmp_mpi' to the target directory
target_dir="$(pwd)/test_simulation"
echo "Copying 'lmp_mpi' to $target_dir"
cp lmp_mpi "$target_dir"

# Check if the copy was successful
if [ $? -ne 0 ]; then
    echo "Error: Copying 'lmp_mpi' failed."
    exit 1
fi

# Step 4: Change directory to 'test_simulation/cac'
cd "$target_dir" || exit

# Step 5: Run 'mpirun'
echo "Running 'mpirun ./lmp_mpi < in_grow.dm'"
mpirun ./lmp_mpi < in_grow.dm

# Check if the 'mpirun' command was successful
if [ $? -ne 0 ]; then
    echo "Error: 'mpirun' command failed."
    exit 1
fi

echo "Build, copy, and simulation completed successfully!"
