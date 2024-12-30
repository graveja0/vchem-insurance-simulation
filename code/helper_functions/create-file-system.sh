#!/usr/bin/env bash
#
# Usage: ./create_directory_structure.sh <project_root_dir>
#
# Example: ./create_directory_structure.sh my_discrete_markov_project
#
# If no argument is provided, it will default to "project_root".

# Get the project root from the command-line argument or default to "project_root"
PROJECT_ROOT="${1:-project_root}"

echo "Creating directory structure in '$PROJECT_ROOT'..."

# Top-level folders
mkdir -p "${PROJECT_ROOT}/data/raw"
mkdir -p "${PROJECT_ROOT}/data/processed"

mkdir -p "${PROJECT_ROOT}/code/0_setup"
mkdir -p "${PROJECT_ROOT}/code/1_data_prep"
mkdir -p "${PROJECT_ROOT}/code/2_multistate_model"
mkdir -p "${PROJECT_ROOT}/code/3_estimation_and_output"
mkdir -p "${PROJECT_ROOT}/code/4_plots_and_viz"
mkdir -p "${PROJECT_ROOT}/code/helper_functions"

mkdir -p "${PROJECT_ROOT}/results/model_objects"
mkdir -p "${PROJECT_ROOT}/results/summary_tables"
mkdir -p "${PROJECT_ROOT}/results/figures"

mkdir -p "${PROJECT_ROOT}/docs/references"
mkdir -p "${PROJECT_ROOT}/docs/methodology"

echo "Directory structure created successfully."